#' Initialization strategy
#'
#' Selecting an initialization point by grouping biomarkers of project and running the SAEM. Initial condition is then selected using maximum log-likelihood.
#'
#' @param project directory of the Monolix project (in .mlxtran). If NULL, the current loaded project is used (default is NULL).
#' @param Nb_genes Size of group of genes.
#' @param unlinkTemporaryProject If temporary project (of genes group) is deleted (defaut: TRUE)
#' @param ObsModel.transfo list containing two lists of transformations and two vectors linking each transformations to their observation model name in the Monolix project. The list should include identity transformations and be named \code{S} and \code{R}. The two vectors should be named \code{linkS} and \code{linkR}.
#'
#' Both \code{S} (for the direct observation models) and \code{linkS}, as well as \code{R} (for latent process models) and \code{linkR}, must have the same length.
#'
#' \itemize{
#'   \item\code{S}: a list of transformations for the direct observation models. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic is observed (from \code{dynFUN});  \item\code{linkS} : a vector specifying the observation model names (that is used in the monolix project, \code{alpha1}, etc.) for each transformation, in the same order as in \code{S};
#'
#'   \item\code{R}: similarly, a list of transformations for the latent process models. Although currently there is only one latent dynamic, each \eqn{s_k, k\leq K} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from \code{dynFUN}; \item \code{linkR} : a vector specifying the observation model names for each transformation, in the same order as in \code{R}.
#' }
#' @param alpha named list of named vector "\code{alpha0}", "\code{alpha1}" (all \code{alpha1} are mandatory). The name of \code{alpha$alpha0} and \code{alpha$alpha1} are the observation model names from the monolix project to which they are linked (if the observations models are defined whithout intercept, alpha$alpha0 need to be set to the vector NULL).
#' @param trueValue -for simulation purposes- named vector of true value for parameters.
#' @param pop.set  population parameters setting for initialization (see \code{\link{setPopulationParameterEstimationSettings}}).
#' @param useSettingsInAPI logical, if the settings for SAEM should not be changed from what is currently set in the project.
#' @param print logical, if the results and algotihm steps should be displayed in the console (default to TRUE).
#' @param digits number of digits to print (default to 2).
#' @param seed value of the seed used to initialize the group (see set.seed).
#' @param conditionalDistributionSampling logical, if conditional distribution estimation should be done on the final project.
#'
#' @returns a list of outputs for every group of genes tested with composition of the group, final parameter estimates, final scores estimates (OFV, AIC, BIC, BICc), temporary project directory. The final selected set is initialize in the project.
#' @export
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' initStrat(project,alpha,ObsModel.transfo,seed=1710)
#' }
initStrat <- function(project,
                      alpha,
                      ObsModel.transfo,
                      Nb_genes = 2,
                      trueValue=NULL,
                      pop.set = NULL,
                      useSettingsInAPI = FALSE,
                      conditionalDistributionSampling = FALSE,
                      print=TRUE,
                      digits=2,
                      unlinkTemporaryProject = TRUE,
                      seed=NULL){

  method <- NULL

  ptm.first <- ptm <- proc.time()
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1
  options(op.new)

  check.proj(project,alpha)
  g = lixoftConnectors::getObservationInformation()
  gy <- data.frame()
  for(var in g$name){
    gy <- rbind(gy,dplyr::rename(g[[var]],"obsid"=var))
  }
  N <- length(unique(gy[["id"]]))
  ntot <- nrow(gy)

  if(useSettingsInAPI){
    pop.set <- lixoftConnectors::getPopulationParameterEstimationSettings()
  }

  pset <- list(nbexploratoryiterations = 200, nbsmoothingiterations = 50,
                simulatedannealing = T, smoothingautostop = T, exploratoryautostop = T)
  if (!is.null(pop.set))
    pset <- modifyList(pset, pop.set[intersect(names(pop.set),
                                                  names(pset))])
  pop.set <- lixoftConnectors::getPopulationParameterEstimationSettings()
  pop.set <- modifyList(pop.set, pset[intersect(names(pset),
                                                   names(pop.set))])

  lixoftConnectors::setPopulationParameterEstimationSettings(pop.set)

  regParam.toprint = paste0(alpha$alpha1,"_pop")
  if(!is.null(trueValue)){
    names(trueValue)[names(trueValue) %in% alpha$alpha1] <- paste0(
      names(trueValue)[names(trueValue) %in% alpha$alpha1],"_pop")
  }

  if(!is.null(alpha$alpha0)){
    para0 = paste0(alpha$alpha0,"_pop")
  }else{
    para0 = c()
  }
  param.toprint = setdiff(union(lixoftConnectors::getPopulationParameterInformation()$name[which(lixoftConnectors::getPopulationParameterInformation()$method=="MLE")],union(para0,sapply(names(alpha$alpha1),FUN=function(yG){lixoftConnectors::getContinuousObservationModel()$parameters[[yG]]},USE.NAMES = FALSE))),regParam.toprint) # the parameter to be estimated minus regParam.toprint

  project.dir <- lixoftConnectors::getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  remix.dir <- file.path(lixoftConnectors::getProjectSettings()$directory,"remix")
  Sys.sleep(0.1)
  if (!dir.exists(remix.dir))
    dir.create(remix.dir)
  Sys.sleep(0.1)
  if (!dir.exists(paste0(remix.dir,"/tmp_init")))
    dir.create(paste0(remix.dir,"/tmp_init"))
  Sys.sleep(0.1)
  summary.file = file.path(remix.dir, "summary_init.txt")
  unlink(summary.file,force=TRUE)
  Sys.sleep(0.1)

  ########################## FIRST ESTIMATION  ###########################
  to.cat <- paste0("\n", dashed.line, "          Starting Initilization by group\n")
  to.cat <- c(to.cat,"   \u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\n")
  print_result(print,summary.file,to.cat,to.print=NULL)

  currentData <- readMLX(project = project,ObsModel.transfo = ObsModel.transfo,alpha = alpha)

  MLE = list()
  for(k in 1:length(alpha$alpha1)){
    yGk = names(alpha$alpha1)[k]
    Yk = currentData$Robs[[yGk]][,yGk]
    if(is.null(alpha$alpha0)){
      MLE[[k]] <- list(estimate=c(mean=0,
                                  sd=sd(Yk)*sqrt((length(Yk)-1)/length(Yk))),
                       sd=c(mean=1e-5,
                            sd=sd(Yk)*sqrt((length(Yk)-1)/length(Yk))/sqrt(2*length(Yk))))
    }else{
      MLE[[k]] <- list(estimate=c(mean=mean(Yk),
                                  sd=sd(Yk)*sqrt((length(Yk)-1)/length(Yk))),
                       sd=c(mean=sd(Yk)*sqrt((length(Yk)-1)/length(Yk))/sqrt(length(Yk)),
                            sd=sd(Yk)*sqrt((length(Yk)-1)/length(Yk))/sqrt(2*length(Yk))))
    }
    names(MLE)[k] <- yGk
  }

  if(!is.null(seed)){
    set.seed(seed)
  }

  Nbr_grp = ceiling(length(alpha$alpha1)/Nb_genes)
  genes = names(alpha$alpha1)
  sizes = rep(length(alpha$alpha1)%/% Nbr_grp,Nbr_grp)
  if((length(alpha$alpha1) %% Nbr_grp)!=0)
    sizes[1:(length(alpha$alpha1)  %% Nbr_grp)] <- sizes[1:(length(alpha$alpha1)  %% Nbr_grp)] + 1
  shuffled = sample(genes)

  split_list <- split(shuffled, rep(1:Nbr_grp, sizes))

  to.cat <- "Group formed :"
  to.cat <- paste0(to.cat,"\n ",paste0(lapply(split_list,FUN=function(x){paste0("(",paste0(sort(x),collapse=","),")")}),collapse=" - "))
  print_result(print,summary.file,to.cat,to.print=NULL)

  res = lapply(split_list,FUN=function(genes){
    to.cat <- ("\n--------------------------------------------------\n")
    to.cat <- paste0(to.cat,"Starting with biomarkers : ",paste0(genes,collapse=", "),"\n")
    print_result(print,summary.file,to.cat,to.print=NULL)

    temporaryDirectory <- paste0(remix.dir,"/tmp_init/project_init_",paste0(genes,collapse="_"),"/")

    if(!dir.exists(temporaryDirectory)){
      dir(temporaryDirectory)
    }

    lixoftConnectors::loadProject(project)
    lixoftConnectors::saveProject(paste0(temporaryDirectory,"project.mlxtran"))

    for(yGk in names(alpha$alpha1)){
      if(!(yGk %in% genes)){
        cmd = paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[[yGk]],"_pop=list(initialValue=0,method='FIXED'))")
        eval(parse(text = cmd))
        if(!is.null(alpha$alpha0)){
          cmd = paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[[yGk]],"_pop=list(initialValue=",MLE[[yGk]]$estimate[["mean"]],",method='FIXED'))")
        eval(parse(text = cmd))
        }
        cmd = paste0("lixoftConnectors::setPopulationParameterInformation(",lixoftConnectors::getContinuousObservationModel()$parameter[[yGk]],"=list(initialValue=",MLE[[yGk]]$estimate[["sd"]],",method='FIXED'))")
        eval(parse(text = cmd))
      }else{
        cmd = paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[[yGk]],"_pop=list(initialValue=1))")
        eval(parse(text = cmd))
      }
    }

    lixoftConnectors::saveProject(paste0(temporaryDirectory,"project.mlxtran"))

    to.cat <- paste0("\nComputing SAEM update for population parameters... \n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    lixoftConnectors::runPopulationParameterEstimation()
    to.cat <- paste0("Estimating log-likelihood... \n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    lixoftConnectors::runLogLikelihoodEstimation()

    lixoftConnectors::saveProject(paste0(temporaryDirectory,"project.mlxtran"))

    to.cat <- "\n > Estimated Population Parameters :\n\n"
    param <- lixoftConnectors::getEstimatedPopulationParameters()

    to.print <- data.frame(EstimatedValue = sapply(param,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})[param.toprint])
    row.names(to.print) <- param.toprint

    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = format(signif(as.numeric(trueValue[param.toprint]),digits=digits),scientific = TRUE),
                        RelativeBias = round((param0[param.toprint]-as.numeric(trueValue[param.toprint]))/as.numeric(trueValue[param.toprint]),digits=digits))
    }
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)

    to.cat <- "\n"
    to.print <- data.frame(EstimatedValue = sapply(param,FUN=function(p){format(signif(p,digits=digits),scientific = T)})[regParam.toprint])
    row.names(to.print) <- regParam.toprint
    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = format(signif(as.numeric(trueValue[regParam.toprint]),digits=digits),scientific = TRUE),
                        RelativeBias = round((param0[regParam.toprint]-as.numeric(trueValue[regParam.toprint]))/as.numeric(trueValue[regParam.toprint]),digits=digits))

      to.print[is.nan(to.print$RelativeBias) | is.infinite(to.print$RelativeBias),"RelativeBias"] <- " "

      to.print[trueValue[regParam.toprint]==0,"TrueValue"] <- "  "
      to.print[param0[regParam.toprint]==0,"EstimatedValue"] <- "  "
    }
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)

    to.cat <- "\n > Estimated logLikelihood :\n"
    to.print <- c(LL=-1/2*lixoftConnectors::getEstimatedLogLikelihood()$importanceSampling[["OFV"]],round(lixoftConnectors::getEstimatedLogLikelihood()$importanceSampling[c("OFV","AIC","BIC","BICc")],digits=digits))
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)


    return(list(genes=genes,
                project=paste0(temporaryDirectory,"project.mlxtran"),
                parameters = param,
                LL=lixoftConnectors::getEstimatedLogLikelihood()$importanceSampling))
  })

  LL = sapply(res,FUN=function(r){r$LL[["OFV"]]})
  keep = which.min(LL)

  to.cat <- ("\n--------------------------------------------------\n")
  to.cat <- paste0(to.cat,"FINAL SET SELECTED : ", paste0(res[[keep]]$genes,collapse=", "),".\n")
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

  lixoftConnectors::loadProject(res[[keep]]$project)
  lixoftConnectors::saveProject(project)

  if(conditionalDistributionSampling){
    to.cat <- "Estimation of the R.E. distribution..."
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    lixoftConnectors::runConditionalDistributionSampling()
  }

  lixoftConnectors::saveProject(project)

  if(unlinkTemporaryProject){
    unlink(paste0(remix.dir,"/tmp_init"),force=TRUE,recursive = TRUE)
  }

  return(res)
}

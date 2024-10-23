#' REMix algorithm
#'
#' @description
#' Regularization and Estimation in MIXed effects model.
#'
#' @details
#' See \code{\link{REMix-package}} for details on the model.
#'
#'
#' @param project directory of the Monolix project (in .mlxtran). If NULL, the current loaded project is used (default is NULL).
#' @param final.project directory of the final Monolix project (default add "_upd" to the Monolix project).
#' @param dynFUN function computing the dynamics of interest for a set of parameters. This function need to contain every sub-function that it may needs (as it is called in a \code{foreach} loop). The output of this function need to return a data.frame with \code{time} as first columns and named dynamics in other columns. It must take in input : \itemize{\item \code{y} a named vector with the initial condition. The names are the dynamics names.
#' \item \code{parms} a named vector of parameter.
#' \item \code{time} vector a timepoint.}
#'
#' See \code{\link{dynFUN_demo}}, \code{\link{model.clairon}}, \code{\link{model.pasin}} or \code{\link{model.pk}} for examples.
#' @param y initial condition of the mechanism model, conform to what is asked in dynFUN.
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
#' @param lambda penalization parameter \eqn{\lambda}.
#' @param eps1 integer (>0) used to define the convergence criteria for the regression parameters.
#' @param eps2 integer (>0) used to define the convergence criteria for the likelihood.
#' @param selfInit logical, if the SAEM is already done in the monolix project should be use as the initial point of the algorithm (if FALSE, SAEM is automatically compute according to \code{pop.set1} settings ; if TRUE, a SAEM through monolix need to have been launched).
#' @param pop.set1 population parameters setting for initialisation (see \code{\link{setPopulationParameterEstimationSettings}}).
#' @param pop.set2 population parameters setting for iterations (see \code{\link{setPopulationParameterEstimationSettings}}).
#' @param pop.set3 population parameters setting for final estimation (see \code{\link{setPopulationParameterEstimationSettings}}).
#' @param prune percentage for prunning (\eqn{\in[0;1]})  in the Adaptative Gauss-Hermite algorithm used to compute the log-likelihood and its derivates (see \code{\link{gh.LL}}).
#' @param n number of points for  gaussian quadrature (see \code{\link{gh.LL}}).
#' @param parallel logical, if the computation should be done in parallel when possible (default TRUE).
#' @param ncores number of cores for parallelization (default NULL and \code{\link{detectCores}} is used).
#' @param print logical, if the results and algotihm steps should be displayed in the console (default to TRUE).
#' @param verbose logical, if progress bar should be printed when possible.
#' @param digits number of digits to print (default to 3).
#' @param trueValue -for simulation purposes- named vector of true value for parameters.
#' @param finalSAEM logical, if a final SAEM should be launch with respect to the final selected set.
#'
#' @seealso \code{\link{cv.remix}}.
#' @return a list of outputs of final project and through the iteration : \itemize{\item \code{info} informations about the parameters (regulatization and population parameter names, alpha names and value of lambda used) ;\item \code{finalRes} containing loglikelihood value \code{LL}, population parameters \code{param} and regularization parameters \code{a.final} values, number of iterations \code{iter} and \code{time} needed, the \code{BIC} of the model built ; \item\code{LL} the vector of Log-Likelihood for the model built over the grid of \eqn{\lambda} ; \item\code{iter.outputs} the list of all remix outputs, i.e. parameters, lieklihood, SAEM estimates and convergence criterion value over the iteration.}
#' @export
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(AB=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,AB=1000)
#' lambda = 1440
#'
#' res = Remix(project = project,
#'             dynFUN = dynFUN_demo,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             eps2=1,
#'             lambda=lambda)
#'
#' plotConvergence(res)
#'
#' trueValue = read.csv(paste0(dirname(project),"/demoSMLX/Simulation/populationParameters.txt"))#'
#'
#' plotSAEM(res,paramToPlot = c("delta_S_pop","phi_S_pop","delta_AB_pop"),trueValue=trueValue)
#' }
remix <- function(project = NULL,
                  final.project = NULL,
                  dynFUN,
                  y,
                  ObsModel.transfo,
                  alpha,
                  lambda,
                  eps1 = 10**(-2),
                  eps2 = 10**(-2),
                  selfInit = FALSE,
                  pop.set1 = NULL,
                  pop.set2 = NULL,
                  pop.set3 = NULL,
                  prune = NULL,
                  n = NULL,
                  parallel=TRUE,
                  ncores = NULL,
                  print = TRUE,
                  verbose=FALSE,
                  digits=3,
                  trueValue = NULL,
                  finalSAEM =FALSE){

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

  ########### Technical PARAMETER ################

  if(parallel){
    if(is.null(ncores)){
      ncores = parallel::detectCores()
    }
    cluster <- snow::makeCluster(ncores)
    doSNOW::registerDoSNOW(cluster)
  }


  ################ START INITIALIZATION OF PROJECT REMIXed ################
  check.proj(project,alpha)

  if(selfInit){
    pop.set1 <- lixoftConnectors::getPopulationParameterEstimationSettings()
  }

  regParam.toprint = paste0(alpha$alpha1,"_pop")
  if(!is.null(trueValue)){
    names(trueValue)[names(trueValue) %in% alpha$alpha1] <- paste0(
      names(trueValue)[names(trueValue) %in% alpha$alpha1],"_pop")
  }
  param.toprint = setdiff(dplyr::filter(lixoftConnectors::getPopulationParameterInformation(),method!="FIXED")$name,regParam.toprint)
  rm.param = setdiff(lixoftConnectors::getPopulationParameterInformation()$name,union(param.toprint,regParam.toprint))

  project.dir <- lixoftConnectors::getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  remix.dir <- file.path(lixoftConnectors::getProjectSettings()$directory,"remix")
  Sys.sleep(0.1)
  if (!dir.exists(remix.dir))
    dir.create(remix.dir)
  Sys.sleep(0.1)
  summary.file = file.path(remix.dir, "summary.txt")
  unlink(summary.file,force=TRUE)
  Sys.sleep(0.1)

  ########################## FIRST ESTIMATION  ###########################
  to.cat <- paste0("\n", dashed.line, " Starting Regulatization and Estimation Algorithm\n")
  to.cat <- c(to.cat,"    \u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\u203E\n")
  print_result(print,summary.file,to.cat,to.print=NULL)

  initial.project <- paste0(remix.dir,"/",sub(".*/([^/]+)/([^/]+)\\.mlxtran", "\\2",project), "_init.mlxtran")

  if (!grepl("\\.", initial.project))
    initial.project <- paste0(initial.project, ".mlxtran")
  if (!grepl("\\.mlxtran", initial.project))
    stop(paste0(initial.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),call. = FALSE)

  if(!selfInit){
    pset1 <- list(nbexploratoryiterations = 200, nbsmoothingiterations = 50,
                  simulatedannealing = F, smoothingautostop = T, exploratoryautostop = T)
    if (!is.null(pop.set1))
      pset1 <- modifyList(pset1, pop.set1[intersect(names(pop.set1),
                                                    names(pset1))])
    pop.set1 <- lixoftConnectors::getPopulationParameterEstimationSettings()
    pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1),
                                                     names(pop.set1))])
  }

  pset2 <- list(nbexploratoryiterations = 150, nbsmoothingiterations = 50, simulatedannealing = F,
                smoothingautostop = F, exploratoryautostop = F)
  if (!is.null(pop.set2))
    pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2),
                                                  names(pset2))])
  pop.set2 <- lixoftConnectors::getPopulationParameterEstimationSettings()
  pop.set2 <- modifyList(pop.set2, pset2[intersect(names(pset2),
                                                   names(pop.set2))])


  suppressMessages({
  check <- check.init(initial.project,pop.set1) # check if initialization step
  if(identical(check,FALSE)){                   # has been done according to
    suppressMessages({                          # the settings given
      lixoftConnectors::loadProject(project)
      lixoftConnectors::saveProject(initial.project)
    })
    check <- check.init(initial.project,pop.set1)
  }
  })

  if(identical(check,FALSE) || !check$SAEM){
    to.cat <- "     - - - Running Initialization Step - - -\n\n"
    to.cat <- c(to.cat,"Initialization using SAEM algorithm :\n")
    to.cat <- c(to.cat,"               -",pop.set1$nbexploratoryiterations,"exploratory iterations;\n")
    to.cat <- c(to.cat,"               -",pop.set1$nbsmoothingiterations,"smoothing iterations.\n\n")
    print_result(print,summary.file,to.cat)
    lixoftConnectors::setPopulationParameterEstimationSettings(pop.set1)
    to.cat <- "Estimation of the population parameters using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    lixoftConnectors::runPopulationParameterEstimation()
    to.cat <- "Estimation of the R.E. distribution using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    lixoftConnectors::runConditionalDistributionSampling()
  }else if(!check$MCMC){
    to.cat <- "     - - - Running Initialization Step - - -\n\n"
    to.cat <- "Estimation of the R.E. distribution using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    lixoftConnectors::runConditionalDistributionSampling()
  }
  lixoftConnectors::saveProject(initial.project)

  if (is.null(final.project)){
    final.project <- paste0(sub(pattern = "(.*)\\..*$",
                                replacement = "\\1", project), "_upd.mlxtran")}
  if (!grepl("\\.", final.project))
    final.project <- paste0(final.project, ".mlxtran")
  if (!grepl("\\.mlxtran", final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),
         call. = FALSE)
  lixoftConnectors::saveProject(final.project)
  if(length(rm.param)==0){
    param0 <- param <- lixoftConnectors::getEstimatedPopulationParameters()
  }else{
    param0 <- param <- lixoftConnectors::getEstimatedPopulationParameters()[-which(names(lixoftConnectors::getEstimatedPopulationParameters())%in%rm.param)]
  }

  ########################## RENDER FIRST ESTIMATION  ###########################
  to.cat <- "\n      - - - <  INITIAL PARAMETERS  > - - -     \n\n"

  to.print <- data.frame(EstimatedValue = sapply(param0,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})[param.toprint])
  row.names(to.print) <- param.toprint
  if(!identical(lixoftConnectors::getEstimatedStandardErrors(),NULL)){
    sd.est = lixoftConnectors::getEstimatedStandardErrors()$stochasticApproximation
    sd.est = sd.est[sd.est$parameter %in% param.toprint,"se"]
    to.print <- cbind(to.print, CI_95 = paste0("[",format(signif(param0[param.toprint]-1.96*sd.est,digits=digits),scientific=TRUE),";",format(signif(param0[param.toprint]+1.96*sd.est,digits=digits),scientific=TRUE),"]"))
  }

  if(!is.null(trueValue)){
    to.print <- cbind(to.print,
                      TrueValue = format(signif(as.numeric(trueValue[param.toprint]),digits=digits),scientific = TRUE),
                      RelativeBias = round((param0[param.toprint]-as.numeric(trueValue[param.toprint]))/as.numeric(trueValue[param.toprint]),digits=digits))
  }
  print_result(print, summary.file, to.cat = to.cat, to.print = to.print)

  to.cat <- "\n"
  to.print <- data.frame(EstimatedValue = sapply(param0,FUN=function(p){format(signif(p,digits=digits),scientific = T)})[regParam.toprint])
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
  ########################## ESTIMATING FIRST LL  ###########################
  to.cat <- "\nEstimating the log-likelihood, and its derivates, using the initial model ... \n"
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  currentData0 <- currentData <-
    readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo,alpha = alpha)
  # FIRST ESTIMATE STORED = NULL

  LL0 <- LL <-
    gh.LL(dynFUN = dynFUN,y = y, data = currentData0, n = n,
          prune = prune, parallel = FALSE)

  LL0$ddLL <- LL$ddLL <- -inflate.H.Ariane(-LL0$ddLL,print=FALSE)
  LL0.pen <- LL.pen <-  LL0$LL - lambda*sum(abs(currentData0$alpha1))


  to.cat <- paste0("             LL : ",round(LL$LL,digits=digits))
  to.cat <- paste0(to.cat,"\n         LL.pen : ",round(LL.pen,digits=digits),"\n")
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

  suppressMessages({lixoftConnectors::loadProject(final.project)})

  a.ini <- a.ini0 <- currentData$alpha1

  ########################## SAVE OUTPUTS   ###########################

  estimates.outputs <- lixoftConnectors::getChartsData("plotSaem")
  LL.outputs <- list(LL0)
  LLpen.outputs <- list(LL0.pen)
  param.outputs <- param0
  crit.outputs <- data.frame()

  stop <- F
  iter =  1
  crit1 <- crit2 <- critb <- 1

  ##########################       ITERATION       ###########################
  while(!stop){
    ############ START ITERATION   ###########
    to.cat <- paste0("   time elapsed : ",round((proc.time()-ptm)["elapsed"],digits=digits),"s\n")
    to.cat <- c(to.cat,dashed.line)
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    to.cat <- c("                 ITERATION ",iter,"\n\n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    ptm <- proc.time()

    ############ UPDATING ALPHA1   ###########
    to.cat <- paste0("Computing taylor update for regularization parameters... \n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    a.final <- setNames(taylorUpdate(alpha = currentData$alpha1,lambda = lambda, dLL = LL0$dLL, ddLL = LL0$ddLL),names(currentData$alpha1))

    currentData$alpha1 <- a.final
    LLpen.aux <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE,onlyLL=TRUE,verbose=verbose) - lambda * sum(abs(a.final))


    if((LLpen.aux %in% c(-Inf,Inf) | LLpen.aux < LL0.pen) && !all(a.final==0)){

      th <- 1e-5
      step <- log(1.5)
      to.recalibrate = which(a.final!=0)
      delta <-  a.final[to.recalibrate] - a.ini[to.recalibrate]

      maxt <- max(abs(delta))

      if(maxt == 0){
        vw <- th
      }else{
        vw <- th/maxt
      }

      res.out.error <- list("old.b" = a.ini[to.recalibrate],
                            "old.rl" = LL0.pen,  # penalise   en l'ancien
                            "old.ca" = critb,    # que les alpha1 ici
                            "old.cb" = crit2)    # pénalise aussi
      # pas tres important, la en cas de plantage, rien d'autres

      sears <- searpas(vw = vw,
                       step = step,
                       b = a.ini[to.recalibrate],
                       delta = delta,
                       funcpa = funcpa,
                       res.out.error = res.out.error,
                       dynFUN=dynFUN,
                       y = y,
                       data = currentData,
                       n = n,
                       prune = prune,
                       stored = NULL,
                       to.recalibrate=to.recalibrate,
                       parallel = FALSE,
                       lambda = lambda)

      a.final <- currentData$alpha1[to.recalibrate]
    }

  to.print <- data.frame(EstimatedValue = format(signif(a.final,digits=digits),scientific=TRUE))
    row.names(to.print) <- regParam.toprint
    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = signif(as.numeric(trueValue[regParam.toprint]),digits=digits),
                        RelativeBias = round(as.numeric((a.final-trueValue[regParam.toprint])/trueValue[regParam.toprint]),digits=digits),
                        " " = ifelse(a.final==0, ifelse(trueValue[regParam.toprint]==0,"  \u2713","  \u2717"),
                                     ifelse(trueValue[regParam.toprint]!=0,"  \u2713","  \u2717"))
      )

      to.print[is.nan(to.print$RelativeBias) | is.infinite(to.print$RelativeBias),"RelativeBias"] <- " "

      to.print[trueValue[regParam.toprint]==0,"TrueValue"] <- "  "
      to.print[a.final==0,"EstimatedValue"] <- "  "
    }
    print_result(print, summary.file, to.cat = NULL, to.print = to.print)


    ############ SAEM UPDATE   ###########
    to.cat <- paste0("\nComputing SAEM update for population parameters... \n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    re <- saemUpdate(project = final.project, final.project = final.project,
                     alpha = alpha, a.final = a.final,iter = iter , pop.set = pop.set2,
                     conditionalDistributionSampling = TRUE, StandardErrors = FALSE)



    estimates = re$SAEMiterations
    for(k in 1:length(alpha$alpha1)){
      cmd = paste0("estimates <- dplyr::mutate(estimates,",alpha$alpha1[k],"_pop =",a.final[k],")")
      eval(parse(text=cmd))
    }

    to.print <- data.frame(EstimatedValue = sapply(re$param,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})[param.toprint])
    row.names(to.print) <- param.toprint
    if(!identical(lixoftConnectors::getEstimatedStandardErrors(),NULL)){
      sd.est = lixoftConnectors::getEstimatedStandardErrors()$stochasticApproximation
      sd.est = sd.est[sd.est$parameter %in% param.toprint,"se"]
      to.print <- cbind(to.print, CI_95 = paste0("[",format(signif(re$param[param.toprint]-1.96*sd.est,digits=digits),scientific=TRUE),";",format(signif(re$param[param.toprint]+1.96*sd.est,digits=digits),scientific=TRUE),"]"))
    }
    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = format(signif(as.numeric(trueValue[param.toprint]),digits=digits),scientific=TRUE),
                        RelativeBias = round(as.numeric((re$param[param.toprint]-trueValue[param.toprint])/trueValue[param.toprint]),digits=digits))
    }
    print_result(print, summary.file, to.cat = NULL, to.print = to.print)

    if(length(rm.param)==0){
      param <- re$param
    }else{
      param <- re$param[-(which(names(re$param) %in% rm.param))]
    }


    ############ ESTIMATE PENALIZED   ###########
    to.cat <- paste0("\nEstimating penalised log-likelihood... \n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    currentData0 <- currentData <- readMLX(project = final.project,
                           ObsModel.transfo = ObsModel.transfo,
                           alpha = alpha)
    LL <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE,verbose=verbose)
    LL$ddLL <- -inflate.H.Ariane(-LL$ddLL,print=FALSE)
    LL.pen <- LL$LL - lambda* sum(abs(a.final))

    to.cat <- paste0("        LL :",round(LL$LL,digits=digits))
    to.cat <- paste0(to.cat,"\n    LL.pen :",round(LL.pen,digits=digits),"\n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    ############ UPDATE CRITERION  ###########
    crit1 = sum(((param0 - param)/ifelse(param0==0,1,param0))**2)
    critb = sum((a.ini0-a.final)**2)
    crit2 = abs(LL0.pen-LL.pen)

    estimates.outputs <- rbind(estimates.outputs,estimates[,colnames(estimates.outputs)])
    LL.outputs <- append(LL.outputs,list(LL))
    LLpen.outputs <- append(LLpen.outputs,LL.pen)
    param.outputs <- rbind(param.outputs,param)

    to.cat <- c("\n\n   Current parameter convergence criterion :",round(crit1,digits=digits),"\n")
    to.cat <- c(to.cat,"  Current  ll pen  convergence criterion  :",round(crit2,digits=digits),"\n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    crit.outputs <- rbind(crit.outputs,data.frame(iter=iter,crit1,critb,crit2))

    if(crit1<eps1 && crit2 <eps2 ){
      stop <- T
    }

    LL0 <- LL
    LL0.pen <- LL.pen
    param0 <- param
    a.ini0 <- a.ini <- a.final
    iter = iter + 1
  }

  if(finalSAEM){
    to.cat <- paste0("   time elapsed : ",round((proc.time()-ptm)["elapsed"],digits=digits),"s\n")
    to.cat <- c(to.cat,dashed.line)
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    to.cat <- c("                 FINAL ITERATION \n\n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    ptm <- proc.time()


    to.cat <- paste0("\nComputing final SAEM... \n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    re <- saemUpdate(project = final.project, final.project = final.project,
                     alpha = alpha, a.final = a.final,iter = iter ,
                     pop.set = pop.set2, pop.setFinal = pop.set3,
                     conditionalDistributionSampling = TRUE,
                     StandardErrors = TRUE, finalSAEM = TRUE )

    ############ ESTIMATE PENALIZED   ###########
    to.cat <- paste0("\nEstimating log-likelihood... \n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    currentData0 <- currentData <- readMLX(project = final.project,
                                           ObsModel.transfo = ObsModel.transfo,
                                           alpha = alpha)
    LLfinal <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE,verbose = verbose,onlyLL = TRUE)

    estimatesfinal = re$SAEMiterations
    for(k in 1:length(alpha$alpha1)){
      if(a.final[k]==0){
        cmd = paste0("estimatesfinal <- dplyr::mutate(estimatesfinal,",alpha$alpha1[k],"_pop =",a.final[k],")")
        eval(parse(text=cmd))
      }
    }

    ################### RENDER FINAL ESTIMATION #########################

    to.cat <- "\n      - - - <  FINAL PARAMETERS  > - - -     \n\n"
    print_result(print,summary.file, to.cat = to.cat,to.print=NULL)


    sd.est = lixoftConnectors::getEstimatedStandardErrors()$stochasticApproximation
    paramtoPrint.FINAL = sd.est$parameter[sd.est$parameter %in% union(regParam.toprint,param.toprint)]
    sd.est = sd.est[sd.est$parameter %in% paramtoPrint.FINAL,"se"]


    to.print <- data.frame(EstimatedValue = sapply(re$param,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})[paramtoPrint.FINAL])
    row.names(to.print) <- paramtoPrint.FINAL
    if(!identical(lixoftConnectors::getEstimatedStandardErrors(),NULL)){
      sd.est = lixoftConnectors::getEstimatedStandardErrors()$stochasticApproximation
      sd.est = sd.est[sd.est$parameter %in% paramtoPrint.FINAL,"se"]
      to.print <- cbind(to.print, CI_95 = paste0("[",format(signif(re$param[paramtoPrint.FINAL]-1.96*sd.est,digits=digits),scientific=TRUE),";",format(signif(re$param[paramtoPrint.FINAL]+1.96*sd.est,digits=digits),scientific=TRUE),"]"))
    }
    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = format(signif(as.numeric(trueValue[paramtoPrint.FINAL]),digits=digits),scientific=TRUE),
                        RelativeBias = round(as.numeric((re$param[paramtoPrint.FINAL]-trueValue[paramtoPrint.FINAL])/trueValue[paramtoPrint.FINAL]),digits=digits))
    }
    print_result(print, summary.file, to.cat = NULL, to.print = to.print)

    if(length(rm.param)==0){
      paramfinal <- re$param
    }else{
      paramfinal <- re$param[-(which(names(re$param) %in% rm.param))]

    }

    to.cat <- "\n - - - <  CRITERION  > - - -     \n"
    to.cat <- paste0(to.cat,"        LL : ",round(LLfinal,digits=digits))
    to.cat <- paste0(to.cat,"\n       BIC :  ",round(-2*LLfinal+log(length(currentData$mu))*sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0),digits=digits),"\n")
    to.cat <- paste0(to.cat,"\n      eBIC :  ",round(-2*LLfinal+log(length(currentData$mu))*sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0)+2*log(choose(length(alpha$alpha1),sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0))),digits=digits),"\n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    ############ outputs  ###########
  }else{
    LLfinal <- LL$LL
    paramfinal <- param

  }

  lixoftConnectors::saveProject(final.project)

  results <- list(info = list(param.toprint=param.toprint,
                              regParam.toprint=regParam.toprint,
                              alpha=alpha,
                              lambda=lambda,
                              finalSAEM = finalSAEM),
                  finalRes=list(LL=LLfinal,
                                param=paramfinal,
                                alpha=paramfinal[paste0(alpha$alpha1,"_pop")],
                                iter=iter,
                                time=(proc.time()-ptm.first)["elapsed"],
                                BIC = -2*LLfinal+log(length(currentData$mu))*sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0),
                                eBIC = -2*LLfinal+log(length(currentData$mu))*sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0)+2*log(choose(length(alpha$alpha1),sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0)))),
                  iterOutputs=list(param=param.outputs,
                                   LL=LL.outputs,
                                   LL.pen = LLpen.outputs,
                                   estimates=estimates.outputs,
                                   criterion = crit.outputs))
  class(results) <- "remix"

  if(parallel)
    snow::stopCluster(cluster)
  return(results)
}


# Inflation matrix from smoothhazard by Bercu -----------------------------
inflate.H.Ariane <- function(H,eps.eigen=10**(-5),print=FALSE){# inflation hessienne ; Merci Ariane ♥

  tr <- sum(diag(H))/ncol(H)

  eigen.values<-eigen(H,symmetric=T,only.values=T)$values

  idpos<-ifelse(any(eigen.values<=eps.eigen),1,0) # nous avons besoin de faire la procédure d'inflation ?
  idpos0<-idpos

  ncount<-0
  if(print){
    cat("BEFORE -----------------\n")
    cat("eigen.values :\n")
    cat("\t",eigen.values,"\n\n")
  }

  while(idpos != 0){
    if(print){cat("\n--------------------------------")
      cat("\nncount =",ncount,"\n")}
    if(ncount==0){ # La c'est le calcul des eta_k dans la formule, faite par Hiviane en pratique dans MLA
      ga <- 0.01
      da <- 1E-2
    }else{
      if(((ncount <= 3) | (ga >= 1)) ){
        da <- da * 5
      }else{
        ga <- ga * 5

        if(ga > 1) ga <- 1
      }
    }
    if(print){
      cat("ga =",ga,"; da =",da,"\n")
    }

    ncount <- ncount + 1

    diagH <- diag(H)
    # La formule de viviane :
    diagH<-ifelse(diagH!=0,diagH+da*(abs((1.e0-ga))*abs(diagH)+ga*tr),
                  da*ga*tr)

    diag(H)<-diagH

    # la on élimine si la essienne a des valeurs infinie (erreur dans ces cas là et on sort de l'algorithme)
    if(sum(H==Inf)>0|sum(H==-Inf)>0){stop("eigen values of hessienne undefined")}
    eigen.values<-eigen(H,symmetric=T,only.values=T)$values
    if(print){
      cat("eigen.values :\n")
      cat("\t",eigen.values,"\n\n")
    }
    # si on a pas réussi, on recommence
    idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
  }

  if(idpos!=0){
    stop("Hessian not defined positive and can't be infalted")
    ite<-ite+1
    pbr_compu<-2
  }
  return(H)}


# Chek --------------------------------------------------------------------
check.init <- function(initial.project,pop.set1){

  initialValue <- NULL

  IndividualParameterModel = lixoftConnectors::getIndividualParameterModel()
  PopulationParameterInformation = lixoftConnectors::getPopulationParameterInformation()
  ContinuousObservationModel = lixoftConnectors::getContinuousObservationModel()

  init.done = tryCatch({prcheck(initial.project)$project
    TRUE },error=function(e){FALSE})
  if(!init.done){
    return(FALSE)
  }else{
    PopulationParameterSettings = lixoftConnectors::getPopulationParameterEstimationSettings()
    if(!identical(PopulationParameterSettings,pop.set1) ||
       !identical(IndividualParameterModel,lixoftConnectors::getIndividualParameterModel()) ||
       !identical(dplyr::select(PopulationParameterInformation,-initialValue),dplyr::select(lixoftConnectors::getPopulationParameterInformation(),-initialValue)) ||
       !identical(ContinuousObservationModel, lixoftConnectors::getContinuousObservationModel())){

      unlink(lixoftConnectors::getProjectSettings()$directory,recursive=TRUE,force=TRUE)
      Sys.sleep(0.1)

      return(FALSE)
    }
  }
  LaunchedTask <- lixoftConnectors::getLaunchedTasks()
  return(list(SAEM=LaunchedTask$populationParameterEstimation,
              MCMC=LaunchedTask$conditionalDistributionSampling))
}

check.proj <- function(project,alpha){
  if (!is.null(project)){
    project <- prcheck(project)$project
  }else{
    project <- lixoftConnectors::getProjectSettings()$project}

  IndividualParameterModel <- lixoftConnectors::getIndividualParameterModel()
  if(any(IndividualParameterModel$distribution[unlist(alpha)]!="normal")){
    cmd = paste0("lixoftConnectors::setIndividualParameterDistribution(",paste0(unlist(alpha),"='normal'",collapse=","),")")
    eval(parse(text=cmd))
    message("[INFO] Distribution of alpha parameters have been switched to normal.")
  }

  ContinuousObservationModel <- lixoftConnectors::getContinuousObservationModel()
  if(any(ContinuousObservationModel$errorModel!="constant")){
    cmd = paste0("lixoftConnectors::setErrorModel(list(",paste0(names(ContinuousObservationModel$errorModel),"='constant'",collapse=","),"))")
    eval(parse(text=cmd))
    message("[INFO] Error Model have been switched to constant.")
  }
  return(invisible(TRUE))
}


# Update ------------------------------------------------------------------
saemUpdate <- function(project = NULL,final.project=NULL,
                       alpha, a.final,
                       iter=NULL,
                       pop.set=NULL,
                       pop.setFinal = NULL,
                       conditionalDistributionSampling = FALSE,
                       StandardErrors = FALSE,
                       finalSAEM = FALSE){

  suppressMessages({
    if (!is.null(project)){
      project <- prcheck(project)$project
    }else{
      project <- lixoftConnectors::getProjectSettings()$project
    }

    lixoftConnectors::loadProject(project)
  })

  if (is.null(final.project)){
    final.project <- paste0(sub(pattern = "(.*)\\..*$",
                                replacement = "\\1", project), "_upd.mlxtran")
  }
  if (!grepl("\\.", final.project))
    final.project <- paste0(final.project, ".mlxtran")
  if (!grepl("\\.mlxtran", final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),
         call. = FALSE)

  if(finalSAEM){
    pset <- list(nbsmoothingiterations=200,nbexploratoryiterations=500,
                 simulatedannealing=T, smoothingautostop=T,exploratoryautostop=T)
    if(!is.null(pop.setFinal))
      pset <-  modifyList(pset, pop.setFinal[intersect(names(pop.setFinal),
                                                       names(pset))])
    pop.set <- lixoftConnectors::getPopulationParameterEstimationSettings()
    pop.set <- modifyList(pop.set, pset[intersect(names(pset), names(pop.set))])
  }else{
    pset <- list(nbsmoothingiterations=50,nbexploratoryiterations=50,
                 simulatedannealing=F, smoothingautostop=F,exploratoryautostop=F)
    if(!is.null(pop.set))
      pset <-  modifyList(pset, pop.set[intersect(names(pop.set),
                                                  names(pset))])
    pop.set <- lixoftConnectors::getPopulationParameterEstimationSettings()
    pop.set <- modifyList(pop.set, pset[intersect(names(pset), names(pop.set))])
  }


  lixoftConnectors::setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)
  if(finalSAEM){
    for(k in 1:length(alpha$alpha1)){
      if(a.final[k]==0){
        eval(parse(text=paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[k],"_pop=list(initialValue=0,method='FIXED'))")))
      }else{
        eval(parse(text=paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[k],"_pop=list(initialValue=",a.final[k],",method='MLE'))")))
      }
    }
  }else{
    for(k in 1:length(alpha$alpha1)){
      eval(parse(text=paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[k],"_pop=list(initialValue=",a.final[k],",method='FIXED'))")))
    }
  }

  lixoftConnectors::setPopulationParameterEstimationSettings(pop.set)

  lixoftConnectors::saveProject(final.project)
  lixoftConnectors::runPopulationParameterEstimation()
  if(conditionalDistributionSampling){
    lixoftConnectors::runConditionalDistributionSampling()
  }
  if(StandardErrors){
    lixoftConnectors::runStandardErrorEstimation()
  }
  lixoftConnectors::saveProject(final.project)

  re = list(SAEMiterations = lixoftConnectors::getChartsData("plotSaem"),
            param = lixoftConnectors::getEstimatedPopulationParameters())
  return(re)
}

taylorUpdate <- function(alpha,lambda,ddLL,dLL){

  # check ; size of alpha and ddLL, dLL

  K = length(alpha)
  X = -ddLL

  value.test = dLL + diag(X)*alpha

  alpha.new = sapply(1:K,FUN = function(k){
    if( value.test[k] < - lambda){
      alpha = alpha[k] + (dLL[k]+lambda)/X[k,k]
    }else if(value.test[k] >  lambda ){
      alpha = alpha[k] + (dLL[k]-lambda)/X[k,k]
    }else{
      alpha = 0
    }
  })

  return(alpha.new)
}

print_result <- function (print, summary.file, to.cat = NULL, to.print = NULL)
{
  if (file.exists(summary.file))
    sink(summary.file, append = TRUE)
  else sink(summary.file)
  if (!is.null(to.cat))
    cat(to.cat)
  if (!is.null(to.print))
    print(to.print)
  sink()
  if (print) {
    if (!is.null(to.cat))
      cat(to.cat)
    if (!is.null(to.print))
      print(to.print)
  }
}

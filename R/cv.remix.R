#' REMixed algorithm over a grid of \eqn{\lambda}
#'
#' @description
#' Regularization and Estimation in MIXed effects model, over a regularization path.
#'
#' @details
#' See \code{\link{REMixed-package}} for details on the model.
#' For each \eqn{\lambda\in\Lambda}, the \code{\link{remix}} is launched.
#' For population parameter estimation settings, see (<https://monolixsuite.slp-software.com/r-functions/2024R1/setpopulationparameterestimationsettings>).
#'
#' @param project directory of the Monolix project (in .mlxtran). If NULL, the current loaded project is used (default is NULL).
#' @param final.project directory of the final Monolix project (default add "_upd" to the Monolix project).
#' @param dynFUN function computing the dynamics of interest for a set of parameters. This function need to contain every sub-function that it may needs (as it is called in a \code{foreach} loop). The output of this function need to return a data.frame with \code{time} as first columns and named dynamics in other columns. It must take in input :
#' \describe{\item{\code{y}}{a named vector with the initial condition. The names are the dynamics names.}
#' \item{\code{parms}}{a named vector of parameter}.
#' \item{\code{time}}{vector a timepoint.}}
#'
#' See \code{\link{dynFUN_demo}}, \code{\link{model.clairon}}, \code{\link{model.pasin}} or \code{\link{model.pk}} for examples.
#' @param y initial condition of the mechanism model, conform to what is asked in dynFUN. If regressor used in Monolix provided a named list of vector of individual initial conditions. Each vector need to be of length 1 (same for all), or exactly the numbre of individuals (range in the same order as their id).
#' @param ObsModel.transfo list containing two lists of transformations and two vectors linking each transformations to their observation model name in the Monolix project. The list should include identity transformations and be named \code{S} and \code{R}. The two vectors should be named \code{linkS} and \code{linkR}.
#'
#' Both \code{S} (for the direct observation models) and \code{linkS}, as well as \code{R} (for latent process models) and \code{linkR}, must have the same length.
#'
#' \describe{
#'   \item{\code{S}}{a list of transformations for the direct observation models. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic is observed (from \code{dynFUN});}\item{\code{linkS}}{a vector specifying the observation model names (that is used in the monolix project, \code{alpha1}, etc.) for each transformation, in the same order as in \code{S}};
#'
#'   \item{\code{R}}{similarly, a list of transformations for the latent process models. Although currently there is only one latent dynamic, each \eqn{s_k, k\leq K} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from \code{dynFUN};} \item{\code{linkR}}{a vector specifying the observation model names for each transformation, in the same order as in \code{R}.}
#' }
#' @param alpha named list of named vector "\code{alpha0}", "\code{alpha1}" (all \code{alpha1} are mandatory). The name of \code{alpha$alpha0} and \code{alpha$alpha1} are the observation model names from the monolix project to which they are linked (if the observations models are defined whithout intercept, alpha$alpha0 need to be set to the vector NULL).
#' @param lambda.grid grid of user-suuplied penalisation parameters for the lasso regularization (if NULL, the sequence is computed based on the data).
#' @param alambda if \code{lambda.grid} is null, coefficients used to compute the grid (default to 0.05, see details).
#' @param nlambda if \code{lambda.grid} is null, number of lambda parameter to test (default to 50).
#' @param lambda_max if \code{lambda.grid} is null, maximum of the lambda grid to test (default is automatically computed, see details)
#' @param eps1 integer (>0) used to define the convergence criteria for the regression parameters.
#' @param eps2 integer (>0) used to define the convergence criteria for the likelihood.
#' @param selfInit logical, if the SAEM is already done in the monolix project should be use as the initial point of the algorithm (if FALSE, SAEM is automatically compute according to \code{pop.set1} settings ; if TRUE, a SAEM through monolix need to have been launched).
#' @param pop.set1 population parameters setting for initialisation (see details).
#' @param pop.set2 population parameters setting for iterations.
#' @param prune percentage for prunning (\eqn{\in[0;1]})  in the Adaptative Gauss-Hermite algorithm used to compute the log-likelihood and its derivates (see \code{\link{gh.LL}}).
#' @param n number of points for  gaussian quadrature (see \code{\link{gh.LL}}).
#' @param parallel logical, if the computation should be done in parallel when possible (default TRUE).
#' @param ncores number of cores for parallelization (default NULL and \code{\link{detectCores}} is used).
#' @param print logical, if the results and algotihm steps should be displayed in the console (default to TRUE).
#' @param digits number of digits to print (default to 3).
#' @param trueValue -for simulation purposes- named vector of true value for parameters.
#' @param unlinkBuildProject logical, if the build project of each lambda should be deleted.
#' @param max.iter maximum number of iteration (default 20).
#'
#' @return
#' A list of outputs of the final project and of the iterative process over each value of \code{lambda.grid}:
#' \describe{
#'   \item{\code{info}}{Information about the parameters.}
#'   \item{\code{project}}{The project path if not unlinked.}
#'   \item{\code{lambda}}{The grid of \eqn{\lambda}.}
#'   \item{\code{BIC}}{Vector of BIC values for the model built over the grid of \eqn{\lambda}.}
#'   \item{\code{BICc}}{Vector of BICc values for the model built over the grid of \eqn{\lambda}.}
#'   \item{\code{LL}}{Vector of log-likelihoods for the model built over the grid of \eqn{\lambda}.}
#'   \item{\code{LL.pen}}{Vector of penalized log-likelihoods for the model built over the grid of \eqn{\lambda}.}
#'   \item{\code{res}}{List of all REMixed results for each \eqn{\lambda} (see \code{\link{remix}}).}
#'   \item{\code{outputs}}{List of all REMixed outputs for each \eqn{\lambda} (see \code{\link{remix}}).}
#' }
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
#'
#' res = cv.remix(project = project,
#'                dynFUN = dynFUN_demo,
#'                y = y,
#'                ObsModel.transfo = ObsModel.transfo,
#'                alpha = alpha,
#'                selfInit = TRUE,
#'                eps1=10**(-2),
#'                ncores=8,
#'                nlambda=8,
#'                eps2=1)
#' }
cv.remix <- function(project = NULL,
                     final.project = NULL,
                     dynFUN,
                     y,
                     ObsModel.transfo,
                     alpha,
                     lambda.grid=NULL,
                     alambda = 0.001,
                     nlambda = 50,
                     lambda_max = NULL,
                     eps1 = 10**(-2),
                     eps2 = 10**(-1),
                     selfInit = FALSE,
                     pop.set1 = NULL,
                     pop.set2 = NULL,
                     prune = NULL,
                     n = NULL,
                     parallel=TRUE,
                     ncores = NULL,
                     print = TRUE,
                     digits=3,
                     trueValue = NULL,
                     unlinkBuildProject = TRUE,
                     max.iter=+Inf){
  method <- NULL

  ptm.first <- ptm <- proc.time()
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  op.original <- options()
  op.new <- options()
  on.exit(options(op.original))
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
  g = mlx.getObservationInformation()
  gy <- data.frame()
  for(var in g$name){
    gy <- rbind(gy,dplyr::rename(g[[var]],"obsid"=var))
  }
  N <- length(unique(gy[["id"]]))
  ntot <- nrow(gy)



  if(selfInit){
    pop.set1 <- mlx.getPopulationParameterEstimationSettings()
  }

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
  param.toprint = setdiff(union(mlx.getPopulationParameterInformation()$name[which(mlx.getPopulationParameterInformation()$method=="MLE")],union(para0,sapply(names(alpha$alpha1),FUN=function(yG){mlx.getContinuousObservationModel()$parameters[[yG]]},USE.NAMES = FALSE))),regParam.toprint) # the parameter to be estimated minus regParam.toprint

  project.dir <- mlx.getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  remix.dir <- file.path(mlx.getProjectSettings()$directory,"remix")
  Sys.sleep(0.1)
  if (!dir.exists(remix.dir))
    dir.create(remix.dir)
  Sys.sleep(0.1)
  summary.file = file.path(remix.dir, "summary_cv.txt")
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
                  simulatedannealing = TRUE, smoothingautostop = TRUE, exploratoryautostop = TRUE)
    if (!is.null(pop.set1))
      pset1 <- modifyList(pset1, pop.set1[intersect(names(pop.set1),
                                                    names(pset1))])
    pop.set1 <- mlx.getPopulationParameterEstimationSettings()
    pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1),
                                                     names(pop.set1))])
  }

  pset2 <- list(nbexploratoryiterations = 150, nbsmoothingiterations = 50, simulatedannealing = TRUE,
                smoothingautostop =TRUE, exploratoryautostop = TRUE)
  if (!is.null(pop.set2))
    pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2),
                                                  names(pset2))])
  pop.set2 <- mlx.getPopulationParameterEstimationSettings()
  pop.set2 <- modifyList(pop.set2, pset2[intersect(names(pset2),
                                                   names(pop.set2))])


  check <- check.init(initial.project,pop.set1) # check if initialization step
  if(identical(check,FALSE)){                   # has been done according to
    suppressMessages({                          # the settings given
      mlx.loadProject(project)
      mlx.saveProject(initial.project)
    })
    check <- check.init(initial.project,pop.set1)
  }

  if(identical(check,FALSE) || !check$SAEM){
    to.cat <- "     - - - Running Initialization Step - - -\n\n"
    to.cat <- c(to.cat,"Initialization using SAEM algorithm :\n")
    to.cat <- c(to.cat,"               -",pop.set1$nbexploratoryiterations,"exploratory iterations;\n")
    to.cat <- c(to.cat,"               -",pop.set1$nbsmoothingiterations,"smoothing iterations.\n\n")
    print_result(print,summary.file,to.cat)
    mlx.setPopulationParameterEstimationSettings(pop.set1)
    to.cat <- "Estimation of the population parameters using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    mlx.runPopulationParameterEstimation()
    to.cat <- "Estimation of the R.E. distribution using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    mlx.runConditionalDistributionSampling()
  }else if(!check$MCMC){
    to.cat <- "     - - - Running Initialization Step - - -\n\n"
    to.cat <- "Estimation of the R.E. distribution using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    mlx.runConditionalDistributionSampling()
  }
  mlx.saveProject(initial.project)

  param0 <- param <- mlx.getEstimatedPopulationParameters()

  ########################## RENDER FIRST ESTIMATION  ###########################
  to.cat <- "\n      - - - <  INITIAL PARAMETERS  > - - -     \n\n"

  to.print <- data.frame(EstimatedValue = sapply(param0,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})[param.toprint])
  row.names(to.print) <- param.toprint

  if(!is.null(trueValue)){
    to.print <- cbind(to.print,
                      TrueValue = format(signif(as.numeric(trueValue[param.toprint]),digits=digits),scientific = TRUE),
                      RelativeBias = round((param0[param.toprint]-as.numeric(trueValue[param.toprint]))/as.numeric(trueValue[param.toprint]),digits=digits))
  }
  print_result(print, summary.file, to.cat = to.cat, to.print = to.print)

  to.cat <- "\n"
  to.print <- data.frame(EstimatedValue = sapply(param0,FUN=function(p){format(signif(p,digits=digits),scientific = TRUE)})[regParam.toprint])
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

  ######################## Computing Regularization PATH   #########################
  to.cat <- "\nComputing regularization path ... \n"
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  currentData0 <- currentData <-
    readMLX(project = initial.project,ObsModel.transfo = ObsModel.transfo,alpha = alpha)


  if(is.null(lambda.grid)){
    if(is.null(lambda_max)){
      lambda_max = lambda.max(dynFUN = dynFUN,y = y, data = currentData0, n = n,
                              prune = prune, parallel=FALSE,verbose=FALSE)
    }
    lambda.grid = lambda_max*(alambda**((1:nlambda)/nlambda))
    lambda.grid <- lambda.grid[lambda.grid!=0]
  }

  # ############### START CV ##########################
  ntasks <- length(lambda.grid)

  keep.lines <- readLines(summary.file)

  cv.res <- lapply(1:length(lambda.grid),FUN=function(array){
    tryCatch({
      prcheck(initial.project)

      lambda = rev(lambda.grid)[array]
      PRINT <- print
      print <- FALSE

      ptm.first <- ptm <- proc.time()

      summary.file.new <- file.path(remix.dir, paste0("summary_",array,".txt"))
      unlink(summary.file.new,force=TRUE)
      Sys.sleep(0.1)
      writeLines(keep.lines,summary.file.new)

      to.cat <- paste0("\n\n Starting algorithm n\u00B0",array,"/",ntasks," with lambda = ",round(lambda,digits=digits),"...\n")
      print_result(PRINT, summary.file, to.cat = to.cat, to.print = NULL)
      print_result(FALSE, summary.file.new, to.cat = to.cat, to.print = NULL)
      to.cat <- paste0("       initialization...\n")
      if(PRINT){cat(to.cat)}

      if (is.null(final.project)){
        final.project <- paste0(remix.dir, paste0("/Build_",array,".mlxtran"))}
      if (!grepl("\\.", final.project))
        final.project <- paste0(final.project, ".mlxtran")
      if (!grepl("\\.mlxtran", final.project))
        stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),
             call. = FALSE)
      mlx.saveProject(final.project)

      ########################## ESTIMATING FIRST LL  ###########################
      to.cat <- "\nEstimating the log-likelihood, and its derivates, using the initial model ... \n"
      print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)

      LL0 <- LL <-
        gh.LL(dynFUN = dynFUN,y = y, data = currentData0, n = n,
              prune = prune, parallel = FALSE,verbose=TRUE)
      LL0$ddLL <- LL$ddLL <- -inflate.H.Ariane(-LL0$ddLL,print=FALSE)
      LL0.pen <- LL.pen <-  LL0$LL - lambda*sum(abs(currentData0$alpha1))


      to.cat <- paste0("             LL : ",round(LL$LL,digits=digits))
      to.cat <- paste0(to.cat,"\n         LL.pen : ",round(LL.pen,digits=digits),"\n")
      print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)

      suppressMessages({mlx.loadProject(final.project)})

      a.ini <- a.ini0 <- currentData$alpha1

      ########################## SAVE OUTPUTS   ###########################

      estimates.outputs <- mlx.getChartsData("plotSaem")
      if(length(setdiff(union(param.toprint,regParam.toprint),colnames(estimates.outputs)[-c(1,2,3)]))!=0){
        for(p in setdiff(union(param.toprint,regParam.toprint),colnames(estimates.outputs)[-c(1,2,3)])){
          estimates.outputs <- cbind(estimates.outputs,new=param[[p]])
          colnames(estimates.outputs)[ncol(estimates.outputs)] <- p
        }
      }
      LL.outputs <- list(LL0)
      LLpen.outputs <- list(LL0.pen)
      param.outputs <- param0
      crit.outputs <- data.frame()

      stop <- FALSE
      iter =  0
      crit1 <- crit2 <- critb <- 1

      ##########################       ITERATION       ###########################
      while(!stop & iter <=max.iter){
        iter <- iter + 1
        to.cat <- paste0("       starting iteration n\u00B0",iter," :\n")
        if(PRINT){cat(to.cat)}
        ############ START ITERATION   ###########
        to.cat <- paste0("   time elapsed : ",round((proc.time()-ptm)["elapsed"],digits=digits),"s\n")
        to.cat <- c(to.cat,dashed.line)
        print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)
        to.cat <- c("                 ITERATION ",iter,"\n\n")
        print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)
        ptm <- proc.time()


        to.cat <- paste0("        update of regulatization parameters...\n")
        if(PRINT){cat(to.cat)}
        ############ UPDATING ALPHA1   ###########
        to.cat <- paste0("Computing taylor update for regularization parameters... \n")
        print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)
        a.final <- setNames(taylorUpdate(alpha = currentData$alpha1,lambda = lambda, dLL = LL0$dLL, ddLL = LL0$ddLL),names(currentData$alpha1))

        currentData$alpha1 <- a.final
        LLpen.aux <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE ,onlyLL=TRUE,verbose = PRINT) - lambda * sum(abs(a.final))

        if((LLpen.aux %in% c(-Inf,Inf) | LLpen.aux < LL0.pen) && !all(a.final==0)){
          to.cat <- "\t/!\ [RECALIBRATION] /!\\n"
          print_result(PRINT, summary.file, to.cat = to.cat, to.print = NULL)
          print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)
          th <- 1
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
                                "old.rl" = LL0.pen,
                                "old.ca" = critb,
                                "old.cb" = crit2)

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
                           print=print,
                           to.recalibrate=to.recalibrate,
                           parallel = FALSE,
                           lambda = lambda)

          a.final[to.recalibrate] <- a.ini[to.recalibrate] + delta*sears$vw

          currentData$alpha1 <- a.final

          LLpen.aux <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE,onlyLL=TRUE,verbose=FALSE) - lambda * sum(abs(a.final))
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
        to.printEND <- to.print
        print_result(print, summary.file.new, to.cat = NULL, to.print = to.print)

        to.cat <- paste0("        update of population parameters...\n")
        if(PRINT){cat(to.cat)}
        ############ SAEM UPDATE   ###########
        to.cat <- paste0("\nComputing SAEM update for population parameters... \n")
        print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)
        re <- saemUpdate(project = final.project, final.project = final.project,
                         currentData = currentData,
                         alpha = alpha, a.final = a.final,iter = iter , pop.set = pop.set2,
                         conditionalDistributionSampling = TRUE, StandardErrors = FALSE)


        estimates = re$SAEMiterations
        if(length(setdiff(union(param.toprint,regParam.toprint),colnames(estimates)[-c(1,2,3)]))!=0){
          for(p in setdiff(union(param.toprint,regParam.toprint),colnames(estimates)[-c(1,2,3)])){
            estimates <- cbind(estimates,new=re$param[[p]])
            colnames(estimates)[ncol(estimates)] <- p
          }
        }

        to.print <- data.frame(EstimatedValue = sapply(re$param,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})[param.toprint])
        row.names(to.print) <- param.toprint
        if(!identical(mlx.getEstimatedStandardErrors(),NULL)){
          sd.est = mlx.getEstimatedStandardErrors()$stochasticApproximation
          sd.est = sd.est[sd.est$parameter %in% param.toprint,"se"]
          to.print <- cbind(to.print, CI_95 = paste0("[",format(signif(re$param[param.toprint]-1.96*sd.est,digits=digits),scientific=TRUE),";",format(signif(re$param[param.toprint]+1.96*sd.est,digits=digits),scientific=TRUE),"]"))
        }
        if(!is.null(trueValue)){
          to.print <- cbind(to.print,
                            TrueValue = format(signif(as.numeric(trueValue[param.toprint]),digits=digits),scientific=TRUE),
                            RelativeBias = round(as.numeric((re$param[param.toprint]-trueValue[param.toprint])/trueValue[param.toprint]),digits=digits))
        }
        print_result(print, summary.file.new, to.cat = NULL, to.print = to.print)
        to.printEND2 <- to.print

        param <- re$param


        ############ ESTIMATE PENALIZED   ###########
        to.cat <- paste0("\nEstimating penalised log-likelihood... \n")
        print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)
        currentData <- readMLX(project = final.project,
                               ObsModel.transfo = ObsModel.transfo,
                               alpha = alpha)
        LL <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE)
        LL$ddLL <- -inflate.H.Ariane(-LL$ddLL,print=FALSE)
        LL.pen <- LL$LL - lambda* sum(abs(a.final))

        to.cat <- paste0("        LL :",round(LL$LL,digits=digits))
        to.cat <- paste0(to.cat,"\n    LL.pen :",round(LL.pen,digits=digits),"\n")
        print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)

        ############ UPDATE CRITERION  ###########
        crit1 = sum(((param0 - param)**2/(param**2+1)))
        critb = sum((a.ini0-a.final)**2)
        crit2 = abs(LL0.pen-LL.pen)

        estimates.outputs <- rbind(estimates.outputs,estimates[,colnames(estimates.outputs)])
        LL.outputs <- append(LL.outputs,list(LL))
        LLpen.outputs <- append(LLpen.outputs,LL.pen)
        param.outputs <- rbind(param.outputs,param)

        to.cat <- c("\n\n   Current parameter convergence criterion :",round(crit1,digits=digits),"\n")
        to.cat <- c(to.cat,"  Current  ll pen  convergence criterion  :",round(crit2,digits=digits),"\n")
        print_result(print, summary.file.new, to.cat = to.cat, to.print = NULL)

        crit.outputs <- rbind(crit.outputs,data.frame(iter=iter,crit1,critb,crit2))

        if(crit1<eps1 && crit2 <eps2 ){
          stop <- TRUE
        }

        LL0 <- LL
        LL0.pen <- LL.pen
        param0 <- param
        a.ini0 <- a.ini <-  a.final
      }
      LLfinal <- LL$LL
      paramfinal <- param

      # mlx.saveProject(final.project)

      results <- list(info = list(project = paste0(remix.dir, paste0("/Build_",array,".mlxtran")),
                                  param.toprint=param.toprint,
                                  regParam.toprint=regParam.toprint,
                                  alpha=alpha,
                                  lambda=lambda,
                                  N=length(currentData$mu),
                                  ntot=ntot),
                      finalRes=list(LL=LLfinal,
                                    LL.pen = LL.pen,
                                    param=paramfinal,
                                    alpha=paramfinal[paste0(alpha$alpha1,"_pop")],
                                    iter=iter,
                                    time=(proc.time()-ptm.first)["elapsed"],
                                    standardError=mlx.getEstimatedStandardErrors()$stochasticApproximation,
                                    saemBeforeTest = NULL),
                      iterOutputs=list(param=param.outputs,
                                       LL=LL.outputs,
                                       LL.pen = LLpen.outputs,
                                       estimates=estimates.outputs,
                                       criterion = crit.outputs))

      Sys.sleep(0.1)
      # progress(array)

      to.cat <- "        DONE !\n"
      to.cat <- "\n      - - - <  CRITERION  > - - -     \n"
      to.cat <- paste0(to.cat,"        LL : ",round(results$finalRes$LL,digits=digits))
      to.cat <- paste0(to.cat,"\n       BIC :  ",round(-2*LLfinal+log(N)*(length(param.toprint)+sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0)),digits=digits))
      to.cat <- paste0(to.cat,"\n      BICc :  ",round(-2*LLfinal+log(ntot)*(sum(paramfinal[paste0(alpha$alpha1,"_pop")]!=0)+PFPR(param.toprint)$PF)+log(N)*PFPR(param.toprint)$PR),"\n")
      print_result(PRINT, summary.file, to.cat = to.cat, to.print = NULL)

      to.cat <- "\n      - - - <   FINAL  PARAMETERS  > - - -     \n\n"
        print_result(PRINT, summary.file, to.cat = to.cat, to.print = to.printEND)
        print_result(PRINT, summary.file, to.cat = NULL, to.print = to.printEND2)
      if(PRINT){cat("\n",dashed.line)}

      return(results)
    },error=function(e){
      message(paste0("Caught an error for lambda =",rev(lambda.grid)[array]," :\n\t>\t ", e$message))
      })
  })

  if(unlinkBuildProject)
    lapply(1:length(lambda.grid),
           FUN=function(array){
             final.project <- paste0(remix.dir, paste0("/Build_",array,".mlxtran"))
             unlink(final.project,recursive=TRUE,force=TRUE)})


  # close(pb)
  if(parallel){
    snow::stopCluster(cluster)
  }

  finalRES = list(info = append(cv.res[[1]]$info[c("param.toprint","regParam.toprint","alpha","N","ntot")],list(project=if(unlinkBuildProject){initial.project}else{sapply(cv.res,FUN=function(f){f$info$project})})),
                  lambda = rev(lambda.grid),
                  LL = sapply(cv.res,FUN=function(f){f$finalRes$LL}),
                  LL.pen = sapply(cv.res,FUN=function(f){f$finalRes$LL - f$info$lambda*sum(abs(as.numeric(f$finalRes$alpha)))}),
                  res = lapply(cv.res,FUN=function(f){f$finalRes}),
                  outputs = lapply(cv.res,FUN=function(f){f$iterOutputs}))

  failed <- which(sapply(finalRES$res,is.null))
  if(length(failed)!=0){
    finalRES$lambda <- finalRES$lambda[-failed]
    # finalRES$BIC <- finalRES$BIC[-failed]
    # finalRES$eBIC <- finalRES$eBIC[-failed]
    finalRES$LL <- finalRES$LL[-failed]
    finalRES$LL.pen <- finalRES$LL.pen[-failed]
    finalRES$res<- finalRES$res[-failed]
    finalRES$outputs<- finalRES$outputs[-failed]
  }

  class(finalRES) <- "cvRemix"
  return(finalRES)
}

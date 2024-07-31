#' REMix algorithm
#'
#' @description
#' Regularization and Estimation in MIXed effects model.
#'
#' @details
#' Suppose that we have a differential system of equations containing variables $(S_{p})_{p\le P}$ and $R$, that depends on some parameters, these dynamics are described by `dynFUN`. We write the process over the time for \eqn{i\le N} individuals, resulting from the differential system for a set of parameters \eqn{\phi_i,(\psi_{li})_{l\le m,i\le N}} for the considered individual \eqn{i\le N}, as \eqn{S_{p}(\cdot,\phi_i,(\psi_{li})_{l\le m})=S_{pi}(\cdot)}, \eqn{p\le P} and \eqn{R(\cdot,\phi_i,(\psi_{li})_{l\le m})=R_i(\cdot)}. Paremeters are described as \deqn{h_l(\psi_{li}) = h_l(\psi_{lpop})+X_i\beta_l + \eta_{li}} with the covariates of individual \eqn{i\le N}, random effects \deqn{\eta_i=(\eta_{li})_{l\le m}\overset{iid}{\sim}\mathcal{N}(\mu_i,\Omega_i)} for \eqn{i\le N} where \eqn{\mu_i} is the estimated random effects of individual \eqn{i} and \eqn{\Omega_i}  is the diagonal matrix of estimated standard deviation of random effects of individual \eqn{i}. The population parameters \eqn{\psi_{pop}=(\psi_{lpop})_{l\le m}}   and \eqn{\beta=(\beta_l)_{l\le m}}  is the vector of covariates effects on parameters.
#' The rest of the population parameters of the structural model, that hasn't random effetcs, are denoted by \eqn{(\phi_i)_{i\le N}}, and are defined as \eqn{\phi_i=\phi_{pop} + X_i \gamma}, they can depends on covariates effects or be constant over the all population.
#' We assume that individual trajectories \eqn{(S_{pi})_{p\le P,i\le N}} are observed through a direct observation model, up to a transformation \eqn{g_p}, \eqn{p\le P}, at differents times \eqn{(t_{pij})_{i\le N,p\le P,j\le n_{ip}}} : \deqn{Y_{pij}=g_p(S_{pi}(t_{pij}))+\epsilon_{pij}} with error \eqn{\epsilon_p=(\epsilon_{pij})\overset{iid}{\sim}\mathcal{N}(0,\varsigma_p^2)} for \eqn{p\le P}.
#' The individual trajectory \eqn{(R_{i})_{i\le N}} is observed through latent processes, up to a transformation \eqn{s_k}, \eqn{k\le K}, observed in \eqn{(t_{kij})_{i\le N,k\le K,j\le n_{kij}}} : \deqn{Z_{kij}=\alpha_{k0}+\alpha_{k1} s_k(R_i(t_{kij}))+\varepsilon_{kij}} where \eqn{\varepsilon_k\overset{iid}{\sim} \mathcal{N}(0,\sigma_k^2)}.
#'
#' @param lambda penalisation parameters for lasso regularization ;
#' @param eps1 convergence limit for parameters distance at each iteraion ;
#' @param eps2 convergence limit for penalized log-likelihood distance at each iteration ;
#' @param pop.set1 Population parameters setting for initialisation ;
#' @param pop.set2 Population parameters setting for algorithm iterations ;
#' @param prune (default NULL) percentage in [0;1] for prunning.
#' @param n (default floor(100**(1/length(theta$psi_pop))) number of points for gaussian quadrature ;
#' @param ncores number of cores for parallelization (default NULL).
#' @param print if TRUE, log are printed in console. Logs are allways saved in a summary file in the remix folder created for the job;
#' @param digits digits to print, (default 3) ;
#' @param trueValue a named vector of true value for parameters (for simulation purpose, if provided, the error is compute at each iteration).
#' @param project the initial Monolix project;
#' @param final.project the final Monolix project (by default adds "_upd" to the original project), every useful log is saved in the remix folder directly in the initial project directory.
#' @param dynFUN Dynamic function ;
#' @param y Initial condition of the model, conform to what is asked in dynFUN ;
#' @param ObsModel.transfo list of 2 list of P,K transformation (need to include identity transformation), named with `S` and `R` : \itemize{\item  ObsModel.transfo$S correspond to the transformation used for direct observation model. For each \eqn{Y_p=h_p(S_p)} the order must be respected and the name indicated which dynamic from dynFUN is observed through this variables \eqn{Y_p}; \item ObsModel.transfo$R correspond to the transformation used for the latent process, as it is now, we only have one latent dynamic so necessarily \eqn{s_k} is applied to `R` but for each \eqn{Y_k} observed, transformation could be different so need to precise as many as in, in same order `alpha$alpha1` ; the name need be set to precise the dynamic from dynFUN to identify the output.}
#' @param alpha named list of named vector "alpha0", "alpha1" (in good order), alpha1 mandatory even if 1. The name of alpha$alpha0 and alpha$alpha1 are the observation model names from the monolix project to which they are linked.
# if alpha_0 vector empty -> all alpha_0 set to 0 ;
# if some alpha_0 not define but not all, set NULL to the one missing
#' @param selfInit If TRUE, the last SAEM done in `project`is used as initialisation for the building algorithm.
#'
#' @return list of outputs of final project and through the iteration
#' @export
#'
#' @seealso \code{\link{cv.Remix}}
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
Remix <- function(project = NULL,
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
                  prune = NULL,
                  n = NULL,
                  ncores = NULL,
                  print = TRUE,
                  digits=3,
                  trueValue = NULL){

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

  if(is.null(ncores)){
    if(.Platform$OS.type == "unix"){
      ncores = parallel::detectCores()
    }else{
      ncores = parallel::detectCores()/2
      if(ncores <1){ncores <- 1}
    }
  }
  cluster <- snow::makeCluster(ncores)
  doSNOW::registerDoSNOW(cluster)


  ################ START INITIALIZATION OF PROJECT REMIXed ################
  check.proj(project,alpha)

  if(selfInit){
    pop.set1 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
  }

  regParam.toprint = paste0(alpha$alpha1,"_pop")
  if(!is.null(trueValue)){
    names(trueValue)[names(trueValue) %in% alpha$alpha1] <- paste0(
      names(trueValue)[names(trueValue) %in% alpha$alpha1],"_pop")
  }
  param.toprint = setdiff(dplyr::filter(Rsmlx:::mlx.getPopulationParameterInformation(),method!="FIXED")$name,regParam.toprint)
  rm.param = setdiff(Rsmlx:::mlx.getPopulationParameterInformation()$name,union(param.toprint,regParam.toprint))

  project.dir <- Rsmlx:::mlx.getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  remix.dir <- file.path(Rsmlx:::mlx.getProjectSettings()$directory,"remix")
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
  Rsmlx:::print_result(print,summary.file,to.cat,to.print=NULL)

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
    pop.set1 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
    pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1),
                                                     names(pop.set1))])
  }

  pset2 <- list(nbexploratoryiterations = 150, nbsmoothingiterations = 50, simulatedannealing = F,
                smoothingautostop = F, exploratoryautostop = F)
  if (!is.null(pop.set2))
    pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2),
                                                  names(pset2))])
  pop.set2 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
  pop.set2 <- modifyList(pop.set2, pset2[intersect(names(pset2),
                                                   names(pop.set2))])


  suppressMessages({
  check <- check.init(initial.project,pop.set1) # check if initialization step
  if(identical(check,FALSE)){                   # has been done according to
    suppressMessages({                          # the settings given
      Rsmlx:::mlx.loadProject(project)
      Rsmlx:::mlx.saveProject(initial.project)
    })
    check <- check.init(initial.project,pop.set1)
  }
  })

  if(identical(check,FALSE) || !check$SAEM){
    to.cat <- "     - - - Running Initialization Step - - -\n\n"
    to.cat <- c(to.cat,"Initialization using SAEM algorithm :\n")
    to.cat <- c(to.cat,"               -",pop.set1$nbexploratoryiterations,"exploratory iterations;\n")
    to.cat <- c(to.cat,"               -",pop.set1$nbsmoothingiterations,"smoothing iterations.\n\n")
    Rsmlx:::print_result(print,summary.file,to.cat)
    Rsmlx:::mlx.setPopulationParameterEstimationSettings(pop.set1)
    to.cat <- "Estimation of the population parameters using the initial model ... \n"
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    Rsmlx:::mlx.runPopulationParameterEstimation()
    to.cat <- "Estimation of the R.E. distribution using the initial model ... \n"
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    Rsmlx:::mlx.runConditionalDistributionSampling()
  }else if(!check$MCMC){
    to.cat <- "     - - - Running Initialization Step - - -\n\n"
    to.cat <- "Estimation of the R.E. distribution using the initial model ... \n"
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    Rsmlx:::mlx.runConditionalDistributionSampling()
  }
  Rsmlx:::mlx.saveProject(initial.project)

  if (is.null(final.project)){
    final.project <- paste0(sub(pattern = "(.*)\\..*$",
                                replacement = "\\1", project), "_upd.mlxtran")}
  if (!grepl("\\.", final.project))
    final.project <- paste0(final.project, ".mlxtran")
  if (!grepl("\\.mlxtran", final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),
         call. = FALSE)
  Rsmlx:::mlx.saveProject(final.project)
  if(length(rm.param)==0){
    param0 <- param <- Rsmlx:::mlx.getEstimatedPopulationParameters()
  }else{
    param0 <- param <- Rsmlx:::mlx.getEstimatedPopulationParameters()[-which(names(Rsmlx:::mlx.getEstimatedPopulationParameters())%in%rm.param)]
  }

  ########################## RENDER FIRST ESTIMATION  ###########################
  to.cat <- "\n      - - - <  INITIAL PARAMETERS  > - - -     \n\n"

  to.print <- data.frame(EstimatedValue = sapply(param0,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})[param.toprint])
  row.names(to.print) <- param.toprint
  if(!identical(Rsmlx:::mlx.getEstimatedStandardErrors(),NULL)){
    sd.est = Rsmlx:::mlx.getEstimatedStandardErrors()$stochasticApproximation
    sd.est = sd.est[sd.est$parameter %in% param.toprint,"se"]
    to.print <- cbind(to.print, CI_95 = paste0("[",format(signif(param0[param.toprint]-1.96*sd.est,digits=digits),scientific=TRUE),";",format(signif(param0[param.toprint]+1.96*sd.est,digits=digits),scientific=TRUE),"]"))
  }
  if(!is.null(trueValue)){
    to.print <- cbind(to.print,
                      TrueValue = format(signif(as.numeric(trueValue[param.toprint]),digits=digits),scientific = TRUE),
                      RelativeBias = round((param0[param.toprint]-as.numeric(trueValue[param.toprint]))/as.numeric(trueValue[param.toprint]),digits=digits))
  }
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)

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
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  ########################## ESTIMATING FIRST LL  ###########################
  to.cat <- "\nEstimating the log-likelihood, and its derivates, using the initial model ... \n"
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
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
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

  suppressMessages({Rsmlx:::mlx.loadProject(final.project)})

  a.ini <- a.ini0 <- currentData$alpha1

  ########################## SAVE OUTPUTS   ###########################

  estimates.outputs <- lixoftConnectors:::getChartsData("plotSaem")
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
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    to.cat <- c("                 ITERATION ",iter,"\n\n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    ptm <- proc.time()

    ############ UPDATING ALPHA1   ###########
    to.cat <- paste0("Computing taylor update for regularization parameters... \n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    a.final <- setNames(taylorUpdate(alpha = currentData$alpha1,lambda = lambda, dLL = LL0$dLL, ddLL = LL0$ddLL),names(currentData$alpha1))

    currentData$alpha1 <- a.final
    LLpen.aux <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE,onlyLL=TRUE) - lambda * sum(abs(a.final))


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
    Rsmlx:::print_result(print, summary.file, to.cat = NULL, to.print = to.print)


    ############ SAEM UPDATE   ###########
    to.cat <- paste0("\nComputing SAEM update for population parameters... \n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
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
    if(!identical(Rsmlx:::mlx.getEstimatedStandardErrors(),NULL)){
      sd.est = Rsmlx:::mlx.getEstimatedStandardErrors()$stochasticApproximation
      sd.est = sd.est[sd.est$parameter %in% param.toprint,"se"]
      to.print <- cbind(to.print, CI_95 = paste0("[",format(signif(re$param[param.toprint]-1.96*sd.est,digits=digits),scientific=TRUE),";",format(signif(re$param[param.toprint]+1.96*sd.est,digits=digits),scientific=TRUE),"]"))
    }
    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = format(signif(as.numeric(trueValue[param.toprint]),digits=digits),scientific=TRUE),
                        RelativeBias = round(as.numeric((re$param[param.toprint]-trueValue[param.toprint])/trueValue[param.toprint]),digits=digits))
    }
    Rsmlx:::print_result(print, summary.file, to.cat = NULL, to.print = to.print)

    param <- re$param[-(which(names(re$param) %in% rm.param))]


    ############ ESTIMATE PENALIZED   ###########
    to.cat <- paste0("\nEstimating penalised log-likelihood... \n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    currentData0 <- currentData <- readMLX(project = final.project,
                           ObsModel.transfo = ObsModel.transfo,
                           alpha = alpha)
    LL <- gh.LL(dynFUN = dynFUN, y = y, data = currentData, n = n, prune = prune, parallel = FALSE)
    LL$ddLL <- -inflate.H.Ariane(-LL$ddLL,print=FALSE)
    LL.pen <- LL$LL - lambda* sum(abs(a.final))

    to.cat <- paste0("        LL :",round(LL$LL,digits=digits))
    to.cat <- paste0(to.cat,"\n    LL.pen :",round(LL.pen,digits=digits),"\n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

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
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

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
  N=length(currentData$mu)

  Rsmlx:::mlx.saveProject(final.project)

  results <- list(info = list(param.toprint=param.toprint,
                              regParam.toprint=regParam.toprint,
                              alpha=alpha,
                              lambda=lambda),
                  finalRes=list(LL=c(Likelihood=LL,PenLikelihood=LL.pen),
                                param=param,
                                alpha=a.final,
                                iter=iter,
                                time=(proc.time()-ptm.first)["elapsed"],
                                BIC = -2*LL0.pen+log(N)*sum(param0[paste0(alpha$alpha1,"_pop")]!=0)),
                  iterOutputs=list(param=param.outputs,
                                   LL=LL.outputs,
                                   LL.pen = LLpen.outputs,
                                   estimates=estimates.outputs,
                                   criterion = crit.outputs))
  class(results) <- "remix"

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
    if(sum(H==Inf)>0|sum(H==-Inf)>0){break}
    eigen.values<-eigen(H,symmetric=T,only.values=T)$values
    if(print){
      cat("eigen.values :\n")
      cat("\t",eigen.values,"\n\n")
    }
    # si on a pas réussi, on recommence
    idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
  }

  if(idpos!=0){
    warning("Hessian not defined positive")
    ite<-ite+1
    pbr_compu<-2
    break
  }
  return(H)}


# Chek --------------------------------------------------------------------
check.init <- function(initial.project,pop.set1){

  IndividualParameterModel = Rsmlx:::mlx.getIndividualParameterModel()
  PopulationParameterInformation = Rsmlx:::mlx.getPopulationParameterInformation()
  ContinuousObservationModel = Rsmlx:::mlx.getContinuousObservationModel()

  init.done = tryCatch({Rsmlx:::prcheck(initial.project)$project
    TRUE },error=function(e){FALSE})
  if(!init.done){
    return(FALSE)
  }else{
    PopulationParameterSettings = Rsmlx:::mlx.getPopulationParameterEstimationSettings()
    if(!identical(PopulationParameterSettings,pop.set1) ||
       !identical(IndividualParameterModel,Rsmlx:::mlx.getIndividualParameterModel()) ||
       !identical(dplyr::select(PopulationParameterInformation,-initialValue),dplyr::select(Rsmlx:::mlx.getPopulationParameterInformation(),-initialValue)) ||
       !identical(ContinuousObservationModel, Rsmlx:::mlx.getContinuousObservationModel())){

      unlink(Rsmlx:::mlx.getProjectSettings()$directory,recursive=TRUE,force=TRUE)
      Sys.sleep(0.1)

      return(FALSE)
    }
  }
  LaunchedTask <- Rsmlx:::mlx.getLaunchedTasks()
  return(list(SAEM=LaunchedTask$populationParameterEstimation,
              MCMC=LaunchedTask$conditionalDistributionSampling))
}

check.proj <- function(project,alpha){
  if (!is.null(project)){
    project <- Rsmlx:::prcheck(project)$project
  }else{
    project <- Rsmlx:::mlx.getProjectSettings()$project}

  IndividualParameterModel <- Rsmlx:::mlx.getIndividualParameterModel()
  if(any(IndividualParameterModel$distribution[unlist(alpha)]!="normal")){
    cmd = paste0("lixoftConnectors::setIndividualParameterDistribution(",paste0(unlist(alpha),"='normal'",collapse=","),")")
    eval(parse(text=cmd))
    message("[INFO] Distribution of alpha parameters have been switched to normal.")
  }

  ContinuousObservationModel <- Rsmlx:::mlx.getContinuousObservationModel()
  if(any(ContinuousObservationModel$errorModel!="constant")){
    cmd = paste0("Rsmlx:::mlx.setErrorModel(list(",paste0(names(ContinuousObservationModel$errorModel),"='constant'",collapse=","),"))")
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
                       conditionalDistributionSampling = FALSE,
                       StandardErrors = FALSE){

  suppressMessages({
    if (!is.null(project)){
      project <- Rsmlx:::prcheck(project)$project
    }else{
      project <- Rsmlx:::mlx.getProjectSettings()$project
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

  pset <- list(nbsmoothingiterations=50,nbexploratoryiterations=50,
               simulatedannealing=F, smoothingautostop=F,exploratoryautostop=F)
  pop.set <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
  pop.set <- modifyList(pop.set, pset[intersect(names(pset), names(pop.set))])

  Rsmlx:::mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)

  for(k in 1:length(alpha$alpha1)){
    eval(parse(text=paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[k],"_pop=list(initialValue=",a.final[k],",method='FIXED'))")))
  }

  Rsmlx:::mlx.setPopulationParameterEstimationSettings(pop.set)

  Rsmlx:::mlx.saveProject(final.project)
  Rsmlx:::mlx.runPopulationParameterEstimation()
  if(conditionalDistributionSampling){
    Rsmlx:::mlx.runConditionalDistributionSampling()
  }
  if(StandardErrors){
    Rsmlx:::mlx.runStandardErrorEstimation()
  }
  Rsmlx:::mlx.saveProject(final.project)

  re = list(SAEMiterations = lixoftConnectors::getChartsData("plotSaem"),
            param = Rsmlx:::mlx.getEstimatedPopulationParameters())
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


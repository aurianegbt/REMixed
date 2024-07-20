#' REMix algorithm
#'
#' @description
#' Regularization and Estimation in MIXed effects model, over a regularization path.
#'
#' @details
#' Suppose that we have a differential system of equations containing variables $(S_{p})_{p\leq P}$ and $R$, that depends on some parameters, these dynamics are described by `dynFun`. We write the process over the time for \eqn{i\leq N} individuals, resulting from the differential system for a set of parameters \eqn{\phi_i,(\psi_{li})_{l\leq m,i\leq N}} for the considered individual \eqn{i\leq N}, as \eqn{S_{p}(\cdot,\phi_i,(\psi_{li})_{l\leq m})=S_{pi}(\cdot)}, \eqn{p\leq P} and \eqn{R(\cdot,\phi_i,(\psi_{li})_{l\leq m})=R_i(\cdot)}. Paremeters are described as \deqn{h_l(\psi_{li}) = h_l(\psi_{lpop})+X_i\beta_l + \eta_{li}} with the covariates of individual \eqn{i\leq N}, random effects \deqn{\eta_i=(\eta_{li})_{l\leq m}\overset{iid}{\sim}\mathcal N(\mu_i,\Omega_i)} for \eqn{i\leq N} where \eqn{\mu_i} is the estimated random effects of individual \eqn{i} and \eqn{\Omega_i}  is the diagonal matrix of estimated standard deviation of random effects of individual \eqn{i}. The population parameters \eqn{\psi_{pop}=(\psi_{lpop})_{l\leq m}}   and \eqn{\beta=(\beta_l)_{l\leq m}}  is the vector of covariates effects on parameters.
#' The rest of the population parameters of the structural model, that hasn't random effetcs, are denoted by \eqn{(\phi_i)_{i\leq N}}, and are defined as \eqn{\phi_i=\phi_{pop} + X_i \gamma}, they can depends on covariates effects or be constant over the all population.
#' We assume that individual trajectories \eqn{(S_{pi})_{p\leq P,i\leq N}} are observed through a direct observation model, up to a transformation \eqn{g_p}, \eqn{p\leq P}, at differents times \eqn{(t_{pij})_{i\leq N,p\leq P,j\leq n_{ip}}} : \deqn{Y_{pij}=g_p(S_{pi}(t_{pij}))+\epsilon_{pij}} with error \eqn{\epsilon_p=(\epsilon_{pij})\overset{iid}{\sim}\mathcal N(0,\varsigma_p^2)} for \eqn{p\leq P}.
#' The individual trajectory \eqn{(R_{i})_{i\leq N}} is observed through latent processes, up to a transformation \eqn{s_k}, \eqn{k\leq K}, observed in \eqn{(t_{kij})_{i\leq N,k\leq K,j\leq n_{kij}}} : \deqn{Z_{kij}=\alpha_{k0}+\alpha_{k1} s_k(R_i(t_{kij}))+\varepsilon_{kij}} where \eqn{\varepsilon_k\overset{iid}{\sim} \mathcal N(0,\sigma_k^2)}.
#'
#' @param project the initial Monolix project;
#' @param final.project the final Monolix project (by default adds "_upd" to the original project), every useful log is saved in the remix folder directly in the initial project directory.
#' @param dynFun Dynamic function ;
#' @param y Initial condition of the model, conform to what is asked in dynFun ;
#' @param ObsModel.transfo list of 2 list of P,K transformation (need to include identity transformation), named with `S` and `R` :
#'
#'   - ObsModel.transfo$S correspond to the transformation used for direct observation model. For each \eqn{Y_p=h_p(S_p)} the order (as in Sobs) must be respected and the name indicated which dynamic from dynFun is observed through this variables \eqn{Y_p};
#'
#'   - ObsModel.transfo$R correspond to the transformation used for the latent process, as it is now, we only have one latent dynamic so necessarily \eqn{s_k} is applied to `R` but for each \eqn{Y_k} observed, transformation could be different so need to precise as many as in Robs ; the name need be set to precise the dynamic from dynFun to identify the output.
#' @param alpha list of named vector "alpha0", "alpha1" (in good order), alpha1 mandatory even if 1
#'
# if alpha_0 vector empty -> all alpha_0 set to 0 ;
# if some alpha_0 not define but not all, set NULL to the one missing
#' @param lambda.grid grid of penalisation parameters for lasso regularization;
#' @param nlambda if lambda.grid is null, number of lambda parameter to use for the lambda.grid;
#' @param eps1 convergence limit for parameters distance at each iteraion ;
#' @param eps2 convergence limit for penalized log-likelihood distance at each iteration ;
#' @param pop.set1 Population parameters setting for initialisation ;
#' @param pop.set2 Population parameters setting for algorithm iterations ;
#' @param prune (default NULL) percentage in [0;1] for prunning.
#' @param n (default floor(100**(1/length(theta$psi_pop))) number of points for gaussian quadrature ;
#' @param ncores number of cores for parallelization (default NULL).
#' @param print if TRUE, log are printed in console. Logs are allways saved in a summary file in the remix folder created for the job;
#' @param digits digits to print, (default 3) ;
#' @param trueValue (FOR SIMULATION, if provided, the error is compute at each iteration. )
#'
#' @return list fo outputs of final project and through the iteration for every lambda on lambda.grid, and the model achieving the best BIC as the best built model.
#' @export
#'
#' @examples
#' [ TO DO ]
cv.Remix <- function(project = NULL,
                  final.project = NULL,
                  dynFUN,
                  y,
                  ObsModel.transfo,
                  alpha,
                  lambda.grid=NULL,
                  nlambda = 50,
                  eps1 = 10**(-1),
                  eps2 = 10**(-1),
                  selfInit = FALSE,
                  pop.set1 = NULL,
                  pop.set2 = NULL,
                  prune = NULL,
                  n = NULL,
                  ncores = NULL,
                  print = TRUE,
                  digits=3,
                  trueValue = NULL){

  ## name0 -> contain the last iterations information
  ## name ->  the current information

  ################ START INITIALIZATION OF PROJECT REMIXed ################

  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1
  options(op.new)

  # load the project
  check.proj(project,alpha) # check if every alpha is normaly distributed,
  # and that each error model is constant

  if(selfInit){
    pop.set1 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
  }

  regParam.toprint = paste0(alpha$alpha1,"_pop")
  if(!is.null(trueValue)){
    names(trueValue)[names(trueValue) %in% alpha$alpha1] <- paste0(
      names(trueValue)[names(trueValue) %in% alpha$alpha1],"_pop")
  }
  param.toprint = setdiff(dplyr::filter(Rsmlx:::mlx.getPopulationParameterInformation(),method!="FIXED")$name,regParam.toprint)
  rm.param = dplyr::filter(Rsmlx:::mlx.getPopulationParameterInformation(),method=="FIXED")$name

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
  to.cat <- c(to.cat,"    ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n")
  Rsmlx:::print_result(print,summary.file,to.cat,to.print=NULL)

  initial.project <- paste0(remix.dir,"/",sub(".*/([^/]+)/([^/]+)\\.mlxtran", "\\2",project), "_init.mlxtran")

  if (!grepl("\\.", initial.project))
    initial.project <- paste0(initial.project, ".mlxtran")
  if (!grepl("\\.mlxtran", initial.project))
    stop(paste0(initial.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),call. = FALSE)

  if(!selfInit){
    pset1 <- list(nbexploratoryiterations = 100, nbsmoothingiterations = 50,
                  simulatedannealing = F, smoothingautostop = T, exploratoryautostop = T)
    if (!is.null(pop.set1))
      pset1 <- modifyList(pset1, pop.set1[intersect(names(pop.set1),
                                                    names(pset1))])
    pop.set1 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
    pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1),
                                                     names(pop.set1))])
  }

  pset2 <- list(nbexploratoryiterations = 50, nbsmoothingiterations = 50, simulatedannealing = F,
                smoothingautostop = F, exploratoryautostop = F)
  if (!is.null(pop.set2))
    pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2),
                                                  names(pset2))])
  pop.set2 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
  pop.set2 <- modifyList(pop.set2, pset2[intersect(names(pset2),
                                                   names(pop.set2))])


  check <- check.init(initial.project,pop.set1) # check if initialization step
  if(identical(check,FALSE)){                   # has been done according to
    suppressMessages({                          # the settings given
      Rsmlx:::mlx.loadProject(project)
      Rsmlx:::mlx.saveProject(initial.project)
    })
    check <- check.init(initial.project,pop.set1)
  }

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

  param0 <- param <- Rsmlx:::mlx.getEstimatedPopulationParameters()[-which(names(Rsmlx:::mlx.getEstimatedPopulationParameters())%in%rm.param)]

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

  ######################## Computing Regularization PATH   #########################
  to.cat <- "\nComputing regularization path ... \n"
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  currentData0 <- currentData <-
    readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo,alpha = alpha)

  if(is.null(lambda.grid)){
    lambda_max = lambda.max(dynFun = dynFun,y = y, data = currentData0, n = n,
                            prune = prune, parallel=FALSE)

    lambda.grid = lambda_max*(0.05**((1:nlambda)/nlambda))
    lambda.grid <- lambda.grid[lambda.grid!=0]
  }




  ############### START CV ##########################
  array = 1
  ntasks <- length(lambda.grid)
  pb <- utils::txtProgressBar(max = ntasks, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cv.res <- foreach::foreach(array = 1:length(lambda.grid),.packages = "REMix",.export = "dynFun",.options.snow=opts)%dopar%{

    Rsmlx:::prcheck(initial.project)


    lambda = lambda.grid[array]
    print = FALSE

    ptm.first <- ptm <- proc.time()

    summary.file.new <- file.path(remix.dir, paste0("summary_",array,".txt"))
    unlink(summary.file.new,force=TRUE)
    Sys.sleep(0.1)
    file.copy(from=summary.file,
              to = summary.file.new)
    summary.file <- summary.file.new

    to.cat <- paste0("array =",array,"\nlambda =",round(lambda.grid[array],digits=digits))
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)


    if (is.null(final.project)){
      final.project <- paste0(remix.dir, paste0("/Build_",array,".mlxtran"))}
    if (!grepl("\\.", final.project))
      final.project <- paste0(final.project, ".mlxtran")
    if (!grepl("\\.mlxtran", final.project))
      stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),
           call. = FALSE)
    Rsmlx:::mlx.saveProject(final.project)

    ########################## ESTIMATING FIRST LL  ###########################
    to.cat <- "\nEstimating the log-likelihood, and its derivates, using the initial model ... \n"
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    LL0 <- LL <-
      gh.LL(dynFun = dynFun,y = y, data = currentData0, n = n,
            prune = prune, parallel=FALSE)
    LL0$ddLL <- LL$ddLL <- -inflate.H.Ariane(-LL0$ddLL,print=FALSE)
    LL0.pen <- LL.pen <-  LL0$LL - lambda*sum(abs(currentData0$alpha1))


    to.cat <- paste0("             LL : ",round(LL$LL,digits=digits))
    to.cat <- paste0(to.cat,"\n         LL.pen : ",round(LL.pen,digits=digits),"\n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    suppressMessages({Rsmlx:::mlx.loadProject(final.project)})

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
      Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
      to.cat <- c("                 ITERATION ",iter,"\n\n")
      Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
      ptm <- proc.time()


      ############ UPDATING ALPHA1   ###########
      to.cat <- paste0("Computing taylor update for regularization parameters... \n")
      Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
      a.final <- taylorUpdate(alpha = currentData$alpha1,lambda = lambda, dLL = LL0$dLL, ddLL = LL0$ddLL)

      currentData$alpha1 <- a.final
      LLpen.aux <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, parallel = FALSE ,onlyLL=TRUE) -lambda * sum(abs(a.final))

      if((LLpen.aux %in% c(-Inf,Inf) | LLpen.aux < LL0.pen) && !all(a.final==0)){
        th <- 1e-5
        step <- log(1.5)
        a.ini[which(a.final==0)] <- 0
        delta <- a.final - a.ini

        maxt <- max(abs(delta))

        if(maxt == 0){
          vw <- th
        }else{
          vw <- th/maxt
        }

        res.out.error <- list("old.b" = a.ini,
                              "old.rl" = LL0.pen,  # pénalisé   en l'ancien
                              "old.ca" = critb,    # que les alpha1 ici
                              "old.cb" = crit2)    # pénalisé aussi

        sears <- marqLevAlg:::searpas(vw = vw,
                                      step = step,
                                      b = a.ini,
                                      delta = delta,
                                      funcpa = funcpa,
                                      res.out.error = res.out.error,
                                      dynFun=dynFun,
                                      y = y,
                                      data = currentData,
                                      n = n,
                                      prune = prune,
                                      parallel=FALSE,
                                      ncores = 1,
                                      onlyLL = TRUE,
                                      lambda = lambda)

        a.final<- a.ini + delta*sears$vw

        currentData$alpha1 <- a.final
      }

      to.print <- data.frame(EstimatedValue = format(signif(a.final,digits=digits),scientific=TRUE))
      row.names(to.print) <- regParam.toprint
      if(!is.null(trueValue)){
        to.print <- cbind(to.print,
                          TrueValue = signif(as.numeric(trueValue[regParam.toprint]),digits=digits),
                          RelativeBias = round(as.numeric((a.final-trueValue[regParam.toprint])/trueValue[regParam.toprint]),digits=digits),
                          " " = ifelse(a.final==0, ifelse(trueValue[regParam.toprint]==0,"  ✓","  ✗"),
                                       ifelse(trueValue[regParam.toprint]!=0,"  ✓","  ✗"))
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
      currentData <- readMLX(project = final.project,
                             ObsModel.transfo = ObsModel.transfo,
                             alpha = alpha)
      LL <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, parallel = FALSE)
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

      to.cat <- c("\n\n   Current parameter criterion :",round(crit1,digits=digits),"\n")
      to.cat <- c(to.cat,"  Current  ll pen  criterion  :",round(crit2,digits=digits),"\n")
      Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

      crit.outputs <- rbind(crit.outputs,data.frame(iter=iter,crit1,critb,crit2))

      if(crit1<eps1 && crit2 <eps2 ){
        stop <- T
      }

      LL0 <- LL
      LL0.pen <- LL.pen
      param0 <- param
      a.ini0 <- a.ini <-  a.final
      iter = iter + 1
      ptm <- proc.time()
    }
    N=length(currentData$mu)

    Rsmlx:::mlx.saveProject(final.project)

    return(list(lambda=lambda,
                res=list(LL=c(Likelihood=LL$LL,PenLikelihood=LL.pen),
                                  param=param,
                                  alpha=a.final,
                                  iter=iter,
                                  time=(proc.time()-ptm.first)["elapsed"],
                                  BIC = -2*LL0.pen+log(N)*sum(param0[paste0(alpha$alpha1,"_pop")]!=0)),
                outputs=list(param=param.outputs,
                             LL=LL.outputs,
                             LL.pen = LLpen.outputs,
                             estimates=estimates.outputs,
                             criterion = crit.outputs)))
  }


  close(pb)
  snow::stopCluster(cluster)
  bestBuild=which.min(sapply(res,FUN=function(r){r$res$res$BIC}))

  return(list(bestBuild=cv.res[[bestBuild]],
              cv.res=cv.res))
}

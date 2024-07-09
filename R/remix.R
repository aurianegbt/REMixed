#' REMix algorithm
#'
#' @description
#' Regularization and Estimation in MIXed effects model.
#'
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
#' @param lambda penalisation parameters for lasso regularization;
#' @param eps1 convergence limit for parameters distance at each iteraion ;
#' @param eps2 convergence limit for penalized log-likelihood distance at each iteration ;
#' @param pop.set1 Population parameters setting for initialisation ;
#' @param pop.set2 Population parameters setting for algorithm iterations ;
#' @param prune (default NULL) percentage in [0;1] for prunning.
#' @param n (default floor(100**(1/length(theta$psi_pop))) number of points for gaussian quadrature ;
#' @param ncores number of cores for parallelization (default NULL).
#' @param print if TRUE, log are printed in console. Logs are allways saved in a summary file in the remix folder created for the job;
#' @param digits digits to print, (default 2) ;
#' @param trueValue (FOR SIMULATION, if provided, the error is compute at each iteration. )
#'
#' @return final project
#' @export
#'
#' @examples
#' # [ TO DO ]
Remix <- function(project = NULL,
                  final.project = NULL,
                  dynFUN,
                  y,
                  ObsModel.transfo,
                  alpha,
                  lambda, eps1 = 10**(-2), eps2 = 10**(-2),
                  selfInit = FALSE,
                  pop.set1 = NULL, pop.set2 = NULL,
                  prune = NULL, n = NULL, ncores = NULL,
                  print = TRUE, digits=2,
                  trueValue = NULL,
                  max.ite.cal = 10){
  # Initialisation projet MLX avec theta^0, alpha^0 déjà faite avant l'appel à cette fonction / on suppose à ce niveau le modèle initialiser, load dans les connecteurs / les modèles d'erreurs doivent tous êtres constants ici....

  # name of observation in monolix y- avec noms de la dyn de dynfun observé, les autres osef mais bien défini dans le bonne ordre !
  # bien mettre le paramètre alpha sur distribution normal
  # bien vérifier modèle d'erreur à constant
  ptm.first <- proc.time()
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1   #hide the warning messages
  options(op.new)

  check.proj(project,alpha)

  if(selfInit){
    pop.set1 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
  }

  param.toprint = dplyr::filter(Rsmlx:::mlx.getPopulationParameterInformation(),method!="FIXED")$name

  project.dir <- Rsmlx:::mlx.getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  remix.dir <- file.path(Rsmlx:::mlx.getProjectSettings()$directory,
                         "remix")
  Sys.sleep(0.1)
  if (!dir.exists(remix.dir))
    dir.create(remix.dir)
  Sys.sleep(0.1)
  summary.file = file.path(remix.dir, "summary.txt")
  unlink(summary.file,force=TRUE)
  Sys.sleep(0.1)

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


  check <- check.init(initial.project,pop.set1)
  if(identical(check,FALSE)){
    suppressMessages({
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
    to.cat <- "Sampling of the conditional distribution using the initial model ... \n"
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    Rsmlx:::mlx.runConditionalDistributionSampling()
  }else if(!check$MCMC){
    to.cat <- "     - - - Running Initialization Step - - -\n\n"
    to.cat <- "Sampling of the conditional distribution using the initial model ... \n"
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
  to.cat <- "\n  -> Estimated parameters :\n"
  to.print <- sapply(Rsmlx:::mlx.getEstimatedPopulationParameters(),FUN=function(p){round(p,digits=5)})[param.toprint]
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  if(!is.null(trueValue)){
    to.cat <- "  -> Error :\n"
    to.print <- sapply(Rsmlx:::mlx.getEstimatedPopulationParameters()-trueValue[names(Rsmlx:::mlx.getEstimatedPopulationParameters())],FUN=function(p){round(p,digits=5)})[param.toprint]
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }


  to.cat <- "\nEstimating the log-likelihood, and its derivates, of the first model ... \n"
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  currentData0 <- currentData <- readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo, alpha = alpha)
  LL0 <- LL <- gh.LL(dynFun = dynFun, y = y, data = currentData0, n = n, prune = prune, ncores = ncores)
  LL0$ddLL <- LL$ddLL <- -inflate.H.Ariane(-LL0$ddLL,print=FALSE)
  LL0.pen <- LL.pen <-  LL0$LL - lambda*sum(abs(currentData0$alpha1))
  to.cat <- paste0("  -> Estimated penalised log-likelihood : ",round(LL.pen,digits=digits),"\n")
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)


  suppressMessages({
    Rsmlx:::mlx.loadProject(final.project)
  })

  estimates.outputs <- lixoftConnectors::getChartsData("plotSaem")
  param0 <- param <- Rsmlx:::mlx.getEstimatedPopulationParameters()

  LL.outputs <- list(LL0)
  LLpen.outputs <- list(LL0$LL-lambda*sum(abs(currentData0$alpha1)))
  param.outputs <- param0
  crit.outputs <- data.frame()

  stop <- F
  iter =  1
  crit1 <- crit2 <- critb <- 1

  while(!stop){
    to.cat <- paste0("  time elapsed :",round((proc.time()-ptm.first)["elapsed"],digits=digits),"\n\n")
    to.cat <- c(to.cat,dashed.line)
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    to.cat <- c("                 ITERATION ",iter,"\n\n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    ptm <- proc.time()
    to.cat <- paste0("Computing taylor update for regularization parameters of ",iter,if(iter==1){"st"}else if(iter==2){"nd"}else if(iter==3){"rd"}else{"th"}," iteration ... \n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    a.ini = currentData$alpha1
    a.final <- taylorUpdate(alpha = currentData$alpha1,lambda = lambda, dLL = LL0$dLL, ddLL = LL0$ddLL)

    # INSERER calibration ici
    currentData$alpha1 <- a.final
    LLpen.aux <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, ncores = ncores,onlyLL=TRUE) -lambda * sum(abs(a.final))

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
                            "old.ca" = critb, # que les beta ici
                            "old.cb" = crit2) # pénalisé aussi

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
                                    ncores = ncores,
                                    onlyLL = TRUE,
                                    lambda = lambda)

      a.new<-a.ini+delta*sears$vw

      currentData$alpha1 <- a.new


      a.final <- a.new
    }

    to.cat <- paste0("Computing SAEM update for population parameters of ",iter,if(iter==1){"st"}else if(iter==2){"nd"}else if(iter==3){"rd"}else{"th"}," iteration ... \n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    re <- saemUpdate(project = final.project, final.project = final.project, alpha = alpha, a.final = a.final,iter = iter , pop.set = pop.set2,conditionalDistributionSampling = TRUE)

    estimates = re$SAEMiterations
    for(k in 1:length(alpha$alpha1)){
      cmd = paste0("estimates <- dplyr::mutate(estimates,",alpha$alpha1[k],"_pop =",a.final[k],")")
      eval(parse(text=cmd))
    }
    to.cat <- "  -> Estimated parameters :\n"
    to.print <- sapply(re$param,FUN=function(p){round(p,digits=5)})[param.toprint]
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
    if(!is.null(trueValue)){
      to.cat <- "\n      Error :\n"
      to.print <- sapply(Rsmlx:::mlx.getEstimatedPopulationParameters()-trueValue[names(Rsmlx:::mlx.getEstimatedPopulationParameters())],FUN=function(p){round(p,digits=5)})[param.toprint]
      Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
    }

    to.cat <- paste0("\nEstimating penalised log-likelihood of ",iter,if(iter==1){"st"}else if(iter==2){"nd"}else if(iter==3){"rd"}else{"th"}," iteration ... \n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    currentData <- readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo, alpha = alpha)
    LL <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, ncores = ncores)
    LL$ddLL <- -inflate.H.Ariane(-LL$ddLL,print=FALSE)
    LL.pen <- LL$LL - lambda* sum(abs(a.final))
    to.cat <- paste0("  -> Estimated penalised log-likelihood :",round(LL.pen,digits=digits),"\n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    param <- re$param

    crit1 = sum((param - param0)**2)
    critb = sum((a.ini-a.final)**2)
    crit2 = abs(LL.pen-LL0.pen)

    estimates.outputs <- rbind(estimates.outputs,estimates[,colnames(estimates.outputs)])
    LL.outputs <- append(LL.outputs,list(LL))
    LLpen.outputs <- append(LLpen.outputs,LL.pen)
    param.outputs <- rbind(param.outputs,param)

    to.cat <- c("\n\n  current parameter criterion:",round(crit1,digits=digits),"\n")
    to.cat <- c(to.cat,"  current penalised ll criterion:",round(crit2,digits=digits),"\n")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)

    crit.outputs <- rbind(crit.outputs,data.frame(iter=iter,crit1,critb,crit2))

    if(crit1<eps1 && crit2 <eps2 ){
      stop <- T
    }

    LL0 <- LL
    LL0.pen <- LL.pen
    param0 <- param
    iter = iter + 1
    ptm <- proc.time()
  }
  N=length(currentData$mu)

  Rsmlx:::mlx.saveProject(final.project)
  return(list(LL0.pen=LL0.pen,
              estimatedParameters=param0,
              BIC=-2*LL0.pen+log(N)*sum(param0[paste0(alpha$alpha1,"_pop")]!=0)))
}

funcpa <- function(b,dynFun, y, data, n, prune, ncores,lambda,onlyLL=TRUE){
  data$alpha1 <- b
  LL = gh.LL(dynFun = dynFun, y = y, data = data, n = n, prune = prune, ncores = ncores,onlyLL=onlyLL)
  return(LL-lambda*sum(abs(b)))
}

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

check.init <- function(initial.project,pop.set1){

  IndividualParameterModel = Rsmlx:::mlx.getIndividualParameterModel()
  PopulationParameterInformation = Rsmlx:::mlx.getPopulationParameterInformation()
  ContinuousObservationModel = Rsmlx:::mlx.getContinuousObservationModel()

  init.done = tryCatch({Rsmlx:::prcheck(initial.project)$project
    TRUE},error=function(e){FALSE})
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
    cmd = paste0("Rsmlx:::mlx.setIndividualParameterDistribution(",paste0(unlist(alpha),"='normal'",collapse=","),")")
    eval(parse(text=cmd))
    message("[INFO] Distribution of alpha parameters have been switched to normal.")
  }

  ContinuousObservationModel <- Rsmlx:::mlx.getContinuousObservationModel()
  if(any(ContinuousObservationModel$errorModel!="constant")){
    cmd = paste0("Rsmlx:::mlx.setErrorModel(",paste0(names(ContinuousObservationModel$errorModel),"='constant'",collapse=","),")")
    eval(parse(text=cmd))
    message("[INFO] Error Model have been switched to constant.")
  }
  return(invisible(TRUE))
}

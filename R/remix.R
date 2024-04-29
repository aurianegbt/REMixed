#' Title
#'
#' @param project
#' @param final.project
#' @param dynFUN
#' @param y
#' @param ObsModel.transfo
#' @param alpha
#' @param lambda
#' @param eps1
#' @param eps2
#' @param pop.set1
#' @param pop.set2
#' @param prune
#' @param n
#' @param ncores
#' @param print
#'
#' @return
#' @export
#'
#' @examples
Remix <- function(project = NULL,
                  final.project = NULL,
                  dynFUN,
                  y,
                  ObsModel.transfo,
                  alpha,
                  lambda, eps1 = 10**(-2), eps2 = 1,
                  pop.set1 = NULL, pop.set2 = NULL,
                  prune = NULL, n = NULL, ncores = NULL,
                  print = TRUE ){
  # Initialisation projet MLX avec theta^0, alpha^0 déjà faite avant l'appel à cette fonction / on suppose à ce niveau le modèle initialiser, load dans les connecteurs / les modèles d'erreurs doivent tous êtres constants ici....

  # name of observation in monolix y- avec noms de la dyn de dynfun observé, les autres osef mais bien défini dans le bonne ordre !
  # bien mettre le paramètre alpha sur distribution normal
  # bien vérifier modèle d'erreur à constant
    ptm <- proc.time()
    dashed.line <- "--------------------------------------------------\n"
    plain.line <- "__________________________________________________\n"
    dashed.short <- "-----------------------\n"
    plain.short <- "_______________________\n"

    op.original <- options()
    op.new <- options()
    op.new$lixoft_notificationOptions$warnings <- 1   #hide the warning messages
    options(op.new)

    check.proj(project,alpha)

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

    to.cat <- paste0("\n", dashed.line, " Starting Regulatization and Estimation Algorithm\n")
    to.cat <- c(to.cat,"    ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n")
    Rsmlx:::print_result(print,summary.file,to.cat,to.print=NULL)

    initial.project <- paste0(remix.dir,"/",sub(".*/([^/]+)/([^/]+)\\.mlxtran", "\\2",project), "_init.mlxtran")

    if (!grepl("\\.", initial.project))
      initial.project <- paste0(initial.project, ".mlxtran")
    if (!grepl("\\.mlxtran", initial.project))
      stop(paste0(initial.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),call. = FALSE)

    pset1 <- list(nbexploratoryiterations = 100, nbsmoothingiterations = 50,
                  simulatedannealing = F, smoothingautostop = T, exploratoryautostop = T)
    if (!is.null(pop.set1))
      pset1 <- modifyList(pset1, pop.set1[intersect(names(pop.set1),
                                                    names(pset1))])
    pop.set1 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
    pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1),
                                                     names(pop.set1))])

    pset2 <- list(nbexploratoryiterations = 50, nbsmoothingiterations = 50, simulatedannealing = F,
                  smoothingautostop = F, exploratoryautostop = F)
    if (!is.null(pop.set2))
      pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2),
                                                    names(pset2))])
    pop.set2 <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
    pop.set2 <- modifyList(pop.set2, pset2[intersect(names(pset2),
                                                     names(pop.set2))])

    check <- check.init(initial.project,pop.set1)
    suppressMessages({
      Rsmlx:::mlx.loadProject(project)
      Rsmlx:::mlx.saveProject(initial.project)
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

    # AVANT GESTION DE PROJET

    currentData0 <- currentData <- readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo, alpha = alpha)
    LL0 <- LL <- gh.LL(dynFun = dynFun, y = y, data = currentData0, n = n, prune = prune, ncores = ncores)
    LL0$ddLL <- LL$ddLL <- -inflate.H.Ariane(-LL0$ddLL)
    LL0.pen <- LL.pen <-  LL0$LL - lambda*sum(abs(currentData0$alpha1))

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

    cat("time")
    print((proc.time()-ptm)["elapsed"])
    ptm <- proc.time()
  while(!stop){
    cat("\n--------------------------------------------------------\n")
    cat("\n start iter",iter,"\n")
    cat("currrent param : \n")
    print(param0)
    cat("current penalised LL :")
    print(LL0.pen)


    a.ini = currentData$alpha1
    a.final <- taylorUpdate(alpha = currentData$alpha1,lambda = lambda, dLL = LL0$dLL, ddLL = LL0$ddLL)

    # INSERER calibration ici
    currentData$alpha1 <- a.final
    LLpen.aux <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, ncores = ncores,onlyLL=TRUE) -lambda * sum(abs(a.final))

    if(LLpen.aux %in% c(-Inf,Inf) | LLpen.aux < LL0.pen){ # test sur la pénalisé !!!!!! pareil pour crit global
      th <- 1e-5
      step <- log(1.5)
      delta <- a.final - a.ini

      maxt <- max(abs(delta))

      if(maxt == 0){
        vw <- th
      }else{
        vw <- th/maxt
      }

      res.out.error <- list("old.b" = a.ini,
                            "old.rl" = LL0.pen,  # pénalisé !!!  en l'ancien
                            "old.ca" = critb, # que les beta ici !!
                            "old.cb" = crit2) # pénalisé aussi !!!!!

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
      LLpen.aux <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, ncores = ncores,onlyLL=TRUE) -lambda * sum(abs(a.final))
    }
    a.final <- a.new

    re <- saemUpdate(project = final.project, final.project = final.project, alpha = alpha, a.final = a.final,iter = iter , pop.set = pop.set2,conditionalDistributionSampling = TRUE)

    estimates = re$SAEMiterations
    for(k in 1:length(alpha$alpha1)){
      cmd = paste0("estimates <- dplyr::mutate(estimates,",alpha$alpha1[k],"_pop =",a.final[k],")")
      eval(parse(text=cmd))
    }

    currentData <- readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo, alpha = alpha)
    LL <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, ncores = ncores)
    LL$ddLL <- -inflate.H.Ariane(-LL$ddLL)
    LL.pen <- LL$LL - sum(abs(a.final))

    param <- re$param

    crit1 = sum((param - param0)**2)
    critb = sum((a.ini-a.final)**2)
    crit2 = abs(LL.pen-LL0.pen)

    estimates.outputs <- rbind(estimates.outputs,estimates[,colnames(estimates.outputs)])
    LL.outputs <- append(LL.outputs,list(LL))
    LLpen.outputs <- append(LLpen.outputs,LL.pen)
    param.outputs <- rbind(param.outputs,param)

    cat("crit param ")
    print(crit1)
    cat("crit LL")
    print(crit2)

    crit.outputs <- rbind(crit.outputs,data.frame(iter=iter,crit1,critb,crit2))

    if(crit1<eps1 && crit2 <eps2 ){
      stop <- T
    }

    LL0 <- LL
    LL0.pen <- LL.pen
    param0 <- param
    iter = iter + 1
    cat("time")
    print((proc.time()-ptm)["elapsed"])
    ptm <- proc.time()
  }

  print(param)
  print(LL)

  # Tâche monolix
  # Lancer un premier SAEM "complet" ??????
  # Lancer Conditional Sampling

  # RENTRER dans l'algo , boucle while

  # Calcul des vraisemblances et dérivées

  # Update alpha

  # Update SAEM

  # conditional sampling monolix pour la suite
  # fin boucle while


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

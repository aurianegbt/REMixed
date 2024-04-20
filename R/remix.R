Remix <- function(project = NULL,
                  final.project = NULL,
                  dynFUN,
                  y,
                  ObsModel.transfo,
                  alpha,
                  lambda, eps1 = 10**(-3), eps2 = 10**(-3),
                  pop.set1 = NULL, pop.set2 = NULL,
                  prune = NULL, n = NULL, ncores = NULL,
                  print = TRUE ){
  # Initialisation projet MLX avec theta^0, alpha^0 déjà faite avant l'appel à cette fonction / on suppose à ce niveau le modèle initialiser, load dans les connecteurs / les modèles d'erreurs doivent tous êtres constants ici....

  # name of observation in monolix y- avec noms de la dyn de dynfun observé, les autres osef mais bien défini dans le bonne ordre !
  # bien mettre le paramètre alpha sur distribution normal
  # bien vérifier modèle d'erreur à constant


  ptm <- proc.time()

  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1   #hide the warning messages
  options(op.new)

  if (!is.null(project)){
    project <- Rsmlx:::prcheck(project)$project
  }else{
    project <- Rsmlx:::mlx.getProjectSettings()$project
  }

  if (is.null(final.project)){
    final.project <- paste0(sub(pattern = "(.*)\\..*$",
                                replacement = "\\1", project), "_upd.mlxtran")}
  if (!grepl("\\.", final.project))
    final.project <- paste0(final.project, ".mlxtran")
  if (!grepl("\\.mlxtran", final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),
         call. = FALSE)

  # CHECK to add

  pset1 <- list(nbexploratoryiterations = 300, nbsmoothingiterations = 200,
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

  # Rsmlx:::mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)

  Rsmlx:::mlx.saveProject(final.project)
  Rsmlx:::mlx.setPopulationParameterEstimationSettings(pop.set1)
  Rsmlx:::mlx.runPopulationParameterEstimation()
  Rsmlx:::mlx.runConditionalDistributionSampling()

  Rsmlx:::mlx.saveProject(final.project)
  currentData0 <- readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo, alpha = alpha)
  LL0 <- gh.LL(dynFun = dynFun, y = y, data = currentData0, n = n, prune = prune, ncores = ncores)

  suppressMessages({
    Rsmlx:::mlx.loadProject(final.project)
  })

  estimates.outputs <- lixoftConnectors::getChartsData("plotSaem")
  param0 = Rsmlx:::mlx.getEstimatedPopulationParameters()

  LL.outputs <- list(LL0)
  param.outputs <- estimates.outputs[nrow(estimates.outputs),-c(1:3),drop=F]

  stop <- F
  iter =  1

  cat("time")
  print((proc.time()-ptm)["elapsed"])
  ptm <- proc.time()
  while(!stop){
    cat("\n start iter",iter,"\n")
    cat("currrent param : \n")
    print(param0)
    cat("current LL :")
    print(LL0$LL)

    a.final <- taylorUpdate(alpha = currentData0$alpha1,lambda = lambda, dLL = LL0$dLL, ddLL = LL0$ddLL)
    re <- saemUpdate(project = final.project, final.project = final.project, alpha = alpha, a.final = a.final,iter = iter , pop.set = pop.set2,conditionalDistributionSampling = TRUE)

    estimates = re$SAEMiterations
    for(k in 1:length(alpha$alpha1)){
      cmd = paste0("estimates <- dplyr::mutate(estimates,",alpha$alpha1[k],"_pop =",a.final[k],")")
      eval(parse(text=cmd))
    }


    currentData <- readMLX(project = final.project,ObsModel.transfo = ObsModel.transfo, alpha = alpha)
    LL <- gh.LL(dynFun = dynFun, y = y, data = currentData, n = n, prune = prune, ncores = ncores)


    estimates.outputs <- rbind(estimates.outputs,estimates[,colnames(estimates.outputs)])
    LL.outputs <- append(LL.outputs,list(LL))
    param <- re$param
    param.outputs <- rbind(param.outputs,param)


    cat("new param : \n")
    print(param)
    cat("new LL :")
    print(LL$LL)


    cat("crit 1")
    print(abs(LL$LL-LL0$LL))
    cat("crit 2 ")
    print(sum((param - param0)**2))

    if(abs(LL$LL-LL0$LL)<eps1 && sum((param - param0)**2)<eps2){
      stop <- T
    }

    LL0 <- LL
    param0 <- param
    iter = iter + 1
    cat("time")
    print((proc.time()-ptm)["elapsed"])
    ptm <- proc.time()
  }

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

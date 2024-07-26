#' Gather information from a monolix project for REMix algorithm.
#'
#' Retrieve all information necessary for the algorithm that can be searched in the created monolix project.
#'
#' @param project directory of the monolix project (in .mlxtran), if NULL get current project load (default to NULL);
#' @param ObsModel.transfo list of 2 list of P,K transformation (it needs to include identity transformation), named with `S` and `R`, and two vector with link to the observation model of the monolix project named with `linkS` and `linkR`. Both `ObsModel.transfo$S` (resp. `ObsModel.transfo$R`) and `ObsModel.transfo$linkS` (resp. `ObsModel.transfo$linkR`) need to have the same length. (see \code{\link{lixoftConnectors::getContinuousObservationModel}} or \code{\link{lixoftConnectors::getObservationInformation}})
#'   - `ObsModel.transfo$S` corresponds to the transformation used for direct observation model. For each \eqn{Y_p=h_p(S_p)}, the name indicate which dynamic from monolix is observed through this variables \eqn{Y_p}. The link vector `ObsModel.transfo$linkS` indicate for each transformation, mention in the same order as in `ObsModel.transfo$S`, what is the observation model name corresponding in monolix.
#'   - `ObsModel.transfo$R` corresponds to the transformation used for the latent process, as it is now, we only have one latent dynamic so necessarily. \eqn{s_k} is applied to the same dynamic but for each \eqn{Y_k} observed, transformation could be different so need to precise as many as in the monolix project created, and the name should be the output from `dynFUN` crresponding. The link vector `ObsModel.transfo$linkR` indicate for each transformation, mention in the same order as in `ObsModel.transfo$R`, what is the observation model name corresponding in monolix.
#' @param alpha list of named vector "alpha0", "alpha1". `alpha$alpha1` is mandatory even if equal to 1. The name of `alpha$alpha0` and `alpha$alpha1` are the observation model names from the monolix project to which they are linked. (if alpha_0 vector empty then all alpha_0 set to 0 ; # if some alpha_0 not define but not all, set NULL to the one missing).
#'
#' @return list with parameters, transformations, observations from the monolix project, with the form needed for the implemented algorithm. If a SAEM hasn't been launched, its the initial condition that are returned.
#'  \itemize{\item `mu` list of estimated mean of the individual random effects distribution a posteriori (if conditional distribution estimation through monolix has been launched for the project).  ;
#'  \item `Omega` list of the estimated variance-covariance matrix of the  individual random effects distribution a posteriori (if conditional distribution estimation through monolix has been launched for the project) ;
#'  \item `theta` list of `phi_pop`, `psi_pop` population parameters (whithout and whith r.e. respectively); `gamma`and `beta` the estimated covariates effects coefficients for each population parameters vector; `omega` the estimated variance of r.e. for `psi_pop` parameters.
#'  \item `alpha1` named vector with the values of regression parameters \eqn{(\alpha_1)_{k\leq K}} linked to our latent process. Names are the observation model names from the monolix project.
#'  \item `covariates`data.frame with covariates (col.) for each individual (row).
#'  \item `ParModel.transfo`and `ParModel.transfo.inv` named list of transformation function for individual parameter model.
#'  \item `Sobs` list (size N) of list of direct observation (size P), each element contain time and observation.
#' \item `Robs` list (size N) of latent observation (size K), each element contain time and observation.
#' \item `Serr` vector of size P containing estimated error model constant.
#' \item `Rerr` vector of size K containing estimated error model constant.
#' \item `ObsModel.transfo` given in input.}
#' @export
#'
#' @seealso \code{\link{Remix}}, \code{\link{cv.Remix}}, \code{\link{lixoftConnectors::getContinuousObservationModel}},  \code{\link{lixoftConnectors::getObservationInformation}}
#'
#' @examples
#' # TODO
#'
readMLX <- function(project=NULL,
                    ObsModel.transfo, # that contains name and info list(S=list(AB=log10),R=list(S=function(x){x},S=function(x){x}))
                    alpha # names,
                    ){

  # check length alpha1 equal to length alpha0 if not null equal to length obs$R
  # check obs output name equal to the one in ObsModel.transfo
  if(length(ObsModel.transfo$S)!=length(ObsModel.transfo$linkS) || length(ObsModel.transfo$R)!=length(ObsModel.transfo$linkR)){
    stop("Error in length of ObsModel.transfo arguments. S and linkS (resp. R and linkR) should have the same length.")
  }
  if(length(alpha$alpha1)!=length(ObsModel.transfo$R)){
    stop("Error in length of alpha1 parameters. alpha$alpha1 and ObsModel.transfo$R should have the same size. ")
  }
  if(!all(names(alpha$alpha1)==ObsModel.transfo$linkR) || (!is.null(alpha$alpha0) && !all(names(alpha$alpha0)==ObsModel.transfo$linkR))){
    stop("Error in names of ObsModel.transfo and alpha. They must be consistent (and in the same order).")
  }

  suppressMessages({
    if (!is.null(project)){
      project <- Rsmlx:::prcheck(project)$project
    }else{
      project <- Rsmlx:::mlx.getProjectSettings()$project
    }

    lixoftConnectors::loadProject(project)
  })

  if(lixoftConnectors::getLaunchedTasks()$populationParameterEstimation){
    est = lixoftConnectors::getEstimatedPopulationParameters()
    value.params = data.frame(name=names(est),initialValue=unname(est))
  }else{
    value.params = lixoftConnectors::getPopulationParameterInformation()[,-3]
    }

  # covariates
  covariates = lixoftConnectors::getCovariateInformation()$covariate[,-1]

  cov.exists = !is.null(covariates)

  # theta
  IndividualParameterModel <- lixoftConnectors::getIndividualParameterModel()
  if(cov.exists){
    formula = IndividualParameterModel$formula
    CovariateModel <- IndividualParameterModel$covariateModel
    if(!ok.beta(formula,CovariateModel)){
      stop("Please take care of having standard covariates coefficients names in your monolix project, of form 'beta_param_cov'.")
    }
    if(!ok.omega(IndividualParameterModel$variability$id,value.params)){
      stop("Please take care of having standard covariates coefficients names in your monolix project, of form 'omega_param'.")
    }
  }


  var.param <- IndividualParameterModel$variability$id
  param <- names(var.param)

  phi.names = param[!var.param & !(param %in% unlist(alpha))]
  psi.names = param[var.param & !(param %in% unlist(alpha))]

  phi_pop = sapply(phi.names,FUN=function(p){
    value.params[value.params$name==paste0(p,"_pop"),"initialValue"]
  })
  psi_pop =  sapply(psi.names,FUN=function(p){
    value.params[value.params$name==paste0(p,"_pop"),"initialValue"]
  })

  if(cov.exists){
    gamma = setNames(lapply(phi.names,FUN=function(p){
      g <-  sapply(names(CovariateModel[[p]]),FUN=function(c){
        if(CovariateModel[[p]][[c]]){
          value.params[value.params$name==paste0("beta_",p,"_",c),
                       "initialValue"]
        }else{
          0
        }},USE.NAMES = F)
      if(identical(g,rep(0,ncol(covariates)))){
        g <- NULL
      }
      return(g)
    }),phi.names)

    beta = setNames(lapply(psi.names,FUN=function(p){
      b <-  sapply(names(CovariateModel[[p]]),FUN=function(c){
        if(CovariateModel[[p]][[c]]){
          value.params[value.params$name==paste0("beta_",p,"_",c),
                       "initialValue"]
        }else{
          0
        }},USE.NAMES = F)
      if(identical(b,rep(0,ncol(covariates)))){
        b <- NULL
      }
      return(b)
    }),psi.names)
  }else{
    gamma = NULL
    beta = NULL
  }

  omega = if(length(psi.names)!=0){
    setNames(sapply(psi.names,FUN=function(p){
      value.params[value.params$name==paste0("omega_",p),"initialValue"]
    }),psi.names)
  }else{NULL}

  alpha0 <- setNames(rep(0,length(ObsModel.transfo$R)),ObsModel.transfo$linkR)
  if(!is.null(alpha$alpha0)){
    alpha0[!is.null(alpha$alpha0)] <- setNames(value.params[value.params$name %in% paste0(alpha$alpha0,"_pop"),"initialValue"],names(alpha$alpha0))
  }
  alpha1 <- setNames(sapply(alpha$alpha1,FUN=function(a){
    value.params[value.params$name==paste0(a,"_pop"),"initialValue"]},USE.NAMES = F),names(alpha$alpha1))

  theta = setNames(list(phi_pop,psi_pop,gamma,beta,alpha0,omega),c("phi_pop","psi_pop","gamma","beta","alpha0","omega"))

  # Serr, Rerr
  ErrorModel <- lixoftConnectors::getContinuousObservationModel()$parameters

  # ObsModel.transfo need to be name after ??? model name in monolix ??? or dyn ??? # need to look at that closer

  Serr <- sapply(ObsModel.transfo$linkS,FUN=function(x){
    value.params[value.params$name==ErrorModel[[x]],"initialValue"]
  },USE.NAMES = T)
  Rerr <- sapply(ObsModel.transfo$linkR,FUN=function(x){
    value.params[value.params$name==ErrorModel[[x]],"initialValue"]
  },USE.NAMES = T)

  # ParModel.transfo and inv
  ParModel <- IndividualParameterModel$distribution[union(phi.names,psi.names)]

  ParModel.transfo = list()
  ParModel.transfo.inv =list()
  for(i in 1:length(ParModel)){
    if(ParModel[[i]]=="normal"){
      ParModel.transfo[[i]] <-  function(x){x}
      ParModel.transfo.inv[[i]] <- function(x){x}
    }else if(ParModel[[i]]=="logNormal" || ParModel[[i]]=="lornormal"){
      ParModel.transfo[[i]] <- log
      ParModel.transfo.inv[[i]] <- exp
    }else if(ParModel[[i]]=="logitNormal"){
      ParModel.transfo[[i]] <- function(x){log(x/(1-x))}
      ParModel.transfo.inv[[i]] <- function(x){exp(x)/(1+exp(x))}
    }
  }
  names(ParModel.transfo) <-  names(ParModel.transfo.inv) <- names(ParModel)




  browser()

  # Sobs Robs
  Observation <- lixoftConnectors::getObservationInformation()

  Sobs <-lapply(ObsModel.transfo$linkS,FUN=function(x){Observation[[x]]})
  Robs <- lapply(ObsModel.transfo$linkR,FUN=function(x){Observation[[x]]})

  if(!lixoftConnectors::getLaunchedTasks()$conditionalDistributionSampling){
    return(list(
      theta = theta,
      alpha1 = alpha1,
      covariates = covariates,
      ParModel.transfo = ParModel.transfo,
      ParModel.transfo.inv = ParModel.transfo.inv,
      Sobs = Sobs,
      Robs = Robs,
      Serr= Serr,
      Rerr=Rerr,
      ObsModel.transfo = ObsModel.transfo
    ))
  }else{
    random.effect = lixoftConnectors::getEstimatedRandomEffects(method="conditionalMean")
    N = length(random.effect$conditionalMean$id)

    Omega = lapply(1:N,FUN=function(i){
      eta = sapply(psi.names,FUN=function(p){
        random.effect$conditionalSD[i,paste0("eta_",p)]
      })
      if(length(eta)==1){
        Omega_i = matrix(eta**2,nrow=1,ncol=1)
      }else{
        Omega_i = diag(eta**2)
      }
      rownames(Omega_i) <- colnames(Omega_i) <- psi.names
      return(Omega_i)
    })

    mu = lapply(1:N,FUN=function(i){
      mu_i = sapply(psi.names,FUN=function(p){
        random.effect$conditionalMean[i,paste0("eta_",p)]
      })
      return(mu_i)
    })
    return(list(
      mu=mu,
      Omega=Omega,
      theta = theta,
      alpha1 = alpha1,
      covariates = covariates,
      ParModel.transfo = ParModel.transfo,
      ParModel.transfo.inv = ParModel.transfo.inv,
      Sobs = Sobs,
      Robs = Robs,
      Serr= Serr,
      Rerr=Rerr,
      ObsModel.transfo = ObsModel.transfo
    ))
  }
}


# OK function -------------------------------------------------------------
ok.beta <- function(formula,CovariateModel){
  form <- strsplit(formula, "\n")[[1]]
  for(p in 1:length(form)){
    f <- form[p]
    cov <- CovariateModel[[p]]
    param <- names(CovariateModel)[p]

    if(!identical(cov,list())){
      for(c in names(cov)[cov]){
        before_c <- sub(paste0("\\*",c,".*"), "", f)
        beta <- sub(".*\\s", "", before_c)
        if(beta != paste0("beta_",param,"_",c)){
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

ok.omega <- function(variability,value.params){
  for(p in names(variability)){
    if(variability[[p]] && !(paste0("omega_",p) %in% value.params$name)){
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Extract Data for REMix Algorithm from Monolix Project
#'
#' This function retrieves all necessary information from a Monolix project file to format the input for the REMix package. It gathers all relevant data required for the REMix algorithm.
#'
#' @param project Directory of the Monolix project (in .mlxtran). If NULL, the current loaded project is used (default is NULL).
#' @param ObsModel.transfo A list containing two lists of transformations and two vectors with their corresponding links to the observation models in the Monolix project. The list should include identity transformations and be named `S` and `R`. It should also include two vectors, `linkS` and `linkR`, which link to the observation models.
#'
#' Both `ObsModel.transfo$S` (for direct observation models) and `ObsModel.transfo$linkS`, as well as `ObsModel.transfo$R` (for latent process models) and `ObsModel.transfo$linkR`, must have the same length.
#'
#' See \code{\link{getContinuousObservationModel}} or \code{\link{getObservationInformation}} for further details.
#'
#' \itemize{
#'   \item `ObsModel.transfo$S`: A list of transformations for the direct observation model. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic from Monolix is observed. The `linkS` vector specifies the observation model name in Monolix for each transformation, in the same order as in `ObsModel.transfo$S`.
#'
#'   \item `ObsModel.transfo$R`: A list of transformations for the latent process model. Although currently there is only one latent dynamic, each \eqn{s_k} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from `dynFUN`. The `linkR` vector specifies the observation model name in Monolix for each transformation, in the same order as in `ObsModel.transfo$R`.
#' }
#' @param alpha A list of named vectors "alpha0" and "alpha1". `alpha$alpha1` is mandatory even if equal to 1. The names in `alpha$alpha0` and `alpha$alpha1` should correspond to the observation model names in the Monolix project. If `alpha$alpha0` is empty, all values are set to 0. If some values in `alpha$alpha0` are missing, they are set to NULL.
#'
#' @return A list containing parameters, transformations, and observations from the Monolix project in the format needed for the REMix algorithm. If an SAEM has not been launched, the initial conditions are returned.
#'
#' \itemize{
#'   \item `mu`: A list of estimated mean values of the individual random effects distribution a posteriori (if conditional distribution estimation through Monolix has been launched).
#'   \item `Omega`: A list of the estimated variance-covariance matrix of the individual random effects distribution a posteriori (if conditional distribution estimation through Monolix has been launched).
#'   \item `theta`: A list of population parameters including `phi_pop`, `psi_pop` (without and with random effects, respectively), `gamma` and `beta` (the estimated covariate effect coefficients for each population parameter vector), and `omega` (the estimated variance of random effects for `psi_pop` parameters).
#'   \item `alpha1`: A named vector with the values of regression parameters \eqn{(\alpha_1)_{k \leq K}} linked to the latent process. Names correspond to the observation model names in the Monolix project.
#'   \item `covariates`: A data.frame with covariates (columns) for each individual (rows).
#'   \item `ParModel.transfo` and `ParModel.transfo.inv`: Named lists of transformation functions for the individual parameter model.
#'   \item `Sobs`: A list (size N) of lists of direct observations (size P), each element containing time and observation data.
#'   \item `Robs`: A list (size N) of latent observations (size K), each element containing time and observation data.
#'   \item `Serr`: A vector of size P containing estimated error model constants.
#'   \item `Rerr`: A vector of size K containing estimated error model constants.
#'   \item `ObsModel.transfo`: The input `ObsModel.transfo` list.
#' }
#' @export
#'
#' @seealso \code{\link{Remix}}, \code{\link{cv.Remix}}, \code{\link{getContinuousObservationModel}}, \code{\link{getObservationInformation}}
#'
#' @examples
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
#' res <- readMLX(project,ObsModel.transfo,alpha)
readMLX <- function(project=NULL,
                    ObsModel.transfo,
                    alpha
                    ){

  check.readMLX(ObsModel.transfo,alpha)

  suppressMessages({
    if (!is.null(project)){
      project <- Rsmlx:::prcheck(project)$project
    }else{
      project <- Rsmlx:::mlx.getProjectSettings()$project
    }

    lixoftConnectors::loadProject(project)
  })

  if(!all(Reduce(union,list(ObsModel.transfo$linkS,ObsModel.transfo$linkR,names(alpha$alpha0),names(alpha$alpha1))) %in% names(lixoftConnectors::getObservationInformation()))){
    stop("Inconsistent observation model between the one provided in the monolix project and the one in function inputs.")
  }

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

  # Sobs Robs
  Observation <- lixoftConnectors::getObservationInformation()

  Sobs <- setNames(lapply(ObsModel.transfo$linkS,FUN=function(x){Observation[[x]]}),ObsModel.transfo$linkS)
  Robs <- setNames(lapply(ObsModel.transfo$linkR,FUN=function(x){Observation[[x]]}),ObsModel.transfo$linkR)

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
check.readMLX <- function(ObsModel.transfo,alpha){
  if(length(ObsModel.transfo$S)!=length(ObsModel.transfo$linkS)){
    stop(paste0("Error in length of ObsModel.transfo arguments. S (length=",length(ObsModel.transfo$S),") and linkS (length=",length(ObsModel.transfo$linkS),") should have the same length."))
  }else if(length(ObsModel.transfo$R)!=length(ObsModel.transfo$linkR)){
    stop(paste0("Error in length of ObsModel.transfo arguments. R (length=",length(ObsModel.transfo$R),") and linkR (length=",length(ObsModel.transfo$linkR),") should have the same length."))
  }else if(length(alpha$alpha1)!=length(ObsModel.transfo$R)){
    stop("Error in length of alpha1 parameters. alpha$alpha1 (length=",length(alpha$alpha1),") and ObsModel.transfo$R (length=",length(ObsModel.transfo$R),") should have the same size. ")
  }else if(!is.null(alpha$alpha0) && length(alpha$alpha1)!=length(ObsModel.transfo$R)){
    stop("Error in length of alpha0 parameters. alpha$alpha0 (length=",length(alpha$alpha0),") and ObsModel.transfo$R (length=",length(ObsModel.transfo$R),") should have the same size. ")
  }else if(!all(names(alpha$alpha1)==ObsModel.transfo$linkR) || (!is.null(alpha$alpha0) && !all(names(alpha$alpha0)==ObsModel.transfo$linkR))){
    stop("Error in names of ObsModel.transfo$linkR and alpha parameters. They must be consistent (and in the same order).")
  }
  return(invisible(TRUE))
}

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

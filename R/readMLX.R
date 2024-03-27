#' Gather information from a monolix project
#'
#' Retrieve all information necessary for the algorithm that can be search in the created monolix project.
#'
#' @param project directory of the monolix project (in .mlxtran), if NULL get current project load, default to NULL
#' @param ObsModel.transfo list of 2 list of P,K transformation (need to include identity transformation), named with `S` and `R` :
#'
#'   - ObsModel.transfo$S correspond to the transformation used for direct observation model. For each \eqn{Y_p=h_p(S_p)}, the name indicate which dynamic from monolix is observed through this variables \eqn{Y_p};
#'
#'   - ObsModel.transfo$R correspond to the transformation used for the latent process, as it is now, we only have one latent dynamic so necessarily \eqn{s_k} is applied to `R` but for each \eqn{Y_k} observed, transformation could be different so need to precise as many as in the monolix project created ;
#' @param alpha list of named vector "alpha0", "alpha1" (in good order), alpha1 mandatory even if 1
#'
# if alpha_0 vector empty -> all alpha_0 set to 0 ;
# if some alpha_0 not define but not all, set NULL to the one missing
#'
#' @return list with parameters, transformations, observations from the monolix project, with the form needed for the implemented algorithm
#' @export
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

  if(is.null(suppressMessages(lixoftConnectors::getLixoftConnectorsState())) || lixoftConnectors::getLixoftConnectorsState()$software!="monolix")
    lixoftConnectors::initializeLixoftConnectors(software = "monolix",force=TRUE)
  if(!is.null(project))
    lixoftConnectors::loadProject(project)

  if(lixoftConnectors::getLaunchedTasks()$populationParameterEstimation){
    est = lixoftConnectors::getEstimatedPopulationParameters()
    value.params = data.frame(name=names(est),initialValue=unname(est))
  }else{
    value.params = lixoftConnectors::getPopulationParameterInformation()[,-3]
    }

  # theta
  IndividualParameterModel <- lixoftConnectors::getIndividualParameterModel()
  formula = IndividualParameterModel$formula
  CovariateModel <- IndividualParameterModel$covariateModel
  if(!ok.beta(formula,CovariateModel)){
    stop("Please take care of having standard covariates coefficients names in your monolix project, of form 'beta_param_cov'.")
  }
  var.param <- IndividualParameterModel$variability$id
  param <- names(var.param)

  phi.names = param[var.param & !(param %in% unlist(alpha))]
  psi.names = param[!var.param & !(param %in% unlist(alpha))]

  phi_pop = sapply(phi.names,FUN=function(p){
    value.params[value.params$name==paste0(p,"_pop"),"initialValue"]
  })
  psi_pop =  sapply(psi.names,FUN=function(p){
    value.params[value.params$name==paste0(p,"_pop"),"initialValue"]
  })

  gamma = setNames(lapply(phi.names,FUN=function(p){
    g <-  sapply(names(CovariateModel[[p]]),FUN=function(c){
      if(CovariateModel[[p]][[c]]){
        value.params[value.params$name==paste0("beta_",p,"_",c),
                     "initialValue"]
      }else{
        0
      }},USE.NAMES = F)
  }),phi.names)

  beta = setNames(lapply(psi.names,FUN=function(p){
    g <-  sapply(names(CovariateModel[[p]]),FUN=function(c){
      if(CovariateModel[[p]][[c]]){
        value.params[value.params$name==paste0("beta_",p,"_",c),
                     "initialValue"]
      }else{
        0
      }},USE.NAMES = F)
  }),psi.names)

  alpha0 <- rep(0,length(ObsModel.transfo$R))
  if(!is.null(alpha$alpha0)){
    alpha0[!is.null(alpha$alpha0)] <- value.params[value.params$name==alpha$alpha0,"initialValue"]
  }
  alpha1 <- sapply(alpha$alpha1,FUN=function(a){
    value.params[value.params$name==paste0(a,"_pop"),"initialValue"]},USE.NAMES = F)

  theta = setNames(list(phi_pop,psi_pop,gamma,beta,alpha0),c("phi_pop","psi_pop","gamma","beta","alpha0"))

  # Serr, Rerr
  ErrorModel <- lixoftConnectors::getContinuousObservationModel()$parameters
  names(ErrorModel) <- stringr::str_remove(names(ErrorModel),"y")

  Serr <- sapply(names(ObsModel.transfo$S),FUN=function(x){
    value.params[value.params$name==ErrorModel[[x]],"initialValue"]
  },USE.NAMES = F)
  Rerr <- sapply(setdiff(names(ErrorModel),names(ObsModel.transfo$S)),FUN=function(x){
    value.params[value.params$name==ErrorModel[[x]],"initialValue"]
  },USE.NAMES = F)

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
  Observation <- lixoftConnectors::getObservationInformation()[paste0("y",names(ErrorModel))]
  names(Observation) <- stringr::str_remove(names(Observation),"y")

  Sobs <-lapply(names(ObsModel.transfo$S),FUN=function(x){Observation[[x]]})
  Robs <- lapply(setdiff(names(Observation),names(ObsModel.transfo$S)),FUN=function(x){Observation[[x]]})

  # covariates
  covariates = lixoftConnectors::getCovariateInformation()$covariate[,-1]

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
    random.effect = getEstimatedRandomEffects(method="conditionalMean")
    N = length(random.effect$conditionalMean$id)

    Omega = lapply(1:N,FUN=function(i){
      eta = sapply(phi.names,FUN=function(p){
        random.effect$conditionalSD[i,paste0("eta_",p)]
      })
      Omega_i = diag(eta**2)
      rownames(Omega_i) <- colnames(Omega_i) <- phi.names
      return(Omega_i)
    })

    mu = lapply(1:N,FUN=function(i){
      mu_i = sapply(phi.names,FUN=function(p){
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

ok.beta <- function(formula,CovariateModel){
  form <- strsplit(formula, "\n")[[1]]
  for(p in 1:length(form)){
    f <- form[p]
    cov <- CovariateModel[[p]]
    param <- names(CovariateModel)[p]

    for(c in names(cov)[cov]){
      before_c <- sub(paste0("\\*",c,".*"), "", f)
      beta <- sub(".*\\s", "", before_c)
      if(beta != paste0("beta_",param,"_",c)){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

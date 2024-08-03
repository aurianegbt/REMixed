#' Extract Data for REMix Algorithm from a Monolix Project
#'
#' @description
#' This function retrieves all necessary information from a Monolix project file to format the input for the REMix package. It gathers all relevant data required for the REMix algorithm.
#'
#' @details
#' To simplify its use, functions \code{\link{remix}}, \code{\link{cv.remix}}, \code{\link{gh.LL}} can be used with arguments \code{data} rather than all necessary informations "\code{theta}", "\code{alpha1}", "\code{covariates}", "\code{ParModel.transfo}", "\code{ParModel.transfo.inv}", "\code{Sobs}", "\code{Robs}", "\code{Serr}", "\code{Rerr}", "\code{ObsModel.transfo}" that could be extract from a monolix project. If the SAEM task of the project hasn't been launched, it's the initial condition and not the estimated parameters that are returned. If the conditional distribution estimation task has been launched, parameters "\code{mu}" and "\code{Omega}" are returned too.
#'
#' @param project directory of the Monolix project (in .mlxtran). If NULL, the current loaded project is used (default is NULL).
#' @param ObsModel.transfo list containing two lists of transformations and two vectors linking each transformations to their observation model name in the Monolix project. The list should include identity transformations and be named \code{S} and \code{R}. The two vectors should be named \code{linkS} and \code{linkR}.
#'
#' Both \code{S} (for the direct observation models) and \code{linkS}, as well as \code{R} (for latent process models) and \code{linkR}, must have the same length.
#'
#' \itemize{
#'   \item\code{S}: a list of transformations for the direct observation models. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic is observed (from \code{dynFUN});  \item\code{linkS} : a vector specifying the observation model names (that is used in the monolix project, \code{alpha1}, etc.) for each transformation, in the same order as in \code{S};
#'
#'   \item\code{R}: similarly, a list of transformations for the latent process models. Although currently there is only one latent dynamic, each \eqn{s_k, k\leq K} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from \code{dynFUN}; \item \code{linkR} : a vector specifying the observation model names for each transformation, in the same order as in \code{R}.
#' }
#' @param alpha named list of named vector "\code{alpha0}", "\code{alpha1}" (all \code{alpha1} are mandatory). The names of \code{alpha$alpha0} and \code{alpha$alpha1} are the observation model names from the monolix project to which they are linked (if the observations models are defined whithout intercept, alpha$alpha0 need to be set to the vector NULL).
#'
#' @return A list containing parameters, transformations, and observations from the Monolix project in the format needed for the REMix algorithm :
#' \itemize{
#'   \item\code{mu} list of individuals random effects estimation (vector of r.e. need to be named by the parameter names), use to locate the density mass (if conditional distribution estimation through Monolix has been launched);
#'   \item\code{Omega} list of individuals estimated standard deviation diagonal matrix (matrix need to have rows and columns named by the parameter names), use to locate the density mass (if conditional distribution estimation through Monolix has been launched);
#'   \item\code{theta} list of model parameters containing i\itemize{\item\code{phi_pop} : named vector with the population parameters with no r.e. \eqn{(\phi_{l\ pop})_{l\leq L}} (NULL if none) ;
#' \item\code{psi_pop} : named vector with the population parameters with r.e. \eqn{(\psi_{l\ pop})_{l\leq m}} ;
#' \item\code{gamma} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with no r.e. ;
#' \item\code{beta} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with r.e..
#' \item\code{alpha0} : named vector of \eqn{(\alpha_{0k})_{k\leq K}} parameters (names are identifier of the observation model, such as in a Monolix project);
#' \item\code{omega} : named vector of estimated r.e. standard deviation;}
#'   \item\code{alpha1}  named vector of regulatization parameters \eqn{(\alpha_{1k})_{k\leq K}}, with identifier of observation model as names;
#'   \item\code{covariates} matrix of individual covariates (size N x n). Individuals must be sorted in the same order than in \code{mu} and \code{Omega};
#'   \item\code{ParModel.transfo} named list of transformation functions \eqn{(h_l)_{l\leq m}} and \eqn{(s_k)_{k\leq K}} for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop}, missing entries are set by default to the identity function ;
#'   \item \code{ParModel.transfo.inv} named list of inverse transformation functions for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop} ;
#'   \item\code{Sobs} ist of individuals trajectories for the direct observation models \eqn{(Y_{pi})_{p \leq P,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{p\leq P} data.frame with time  \eqn{(t_{pij})_{j\leq n_{ip}}} and observations  \eqn{(Y_{pij})_{j\leq n_{ip}}}. Each data.frame is named with the observation model identifiers ;
#'   \item\code{Robs} list of individuals trajectories for the latent observation models \eqn{(Z_{ki})_{k \leq K,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{k\leq K} data.frame with time  \eqn{(t_{kij})_{j\leq n_{ik}}} and observations  \eqn{(Z_{kij})_{j\leq n_{ik}}}. Each data.frame is named with the observation model identifiers ;
#'   \item\code{Serr}  named vector of the estimated error mocel constants \eqn{(\varsigma_p)_{p\leq P}} with observation model identifiers as names ;
#'   \item\code{Rerr} named vector of the estimated error mocel constants \eqn{(\sigma_k)_{k\leq K}} with observation model identifiers as names ;
#'   \item\code{ObsModel.transfo} same as input\code{ObsModel.transfo} list.
#' }
#' @export
#'
#' @seealso \code{\link{remix}}, \code{\link{cv.remix}}, \code{\link{getContinuousObservationModel}}, \code{\link{getObservationInformation}}
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

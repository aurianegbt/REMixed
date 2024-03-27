#' Generate individual parameters
#'
#' @description
#' Generate the individual parameters of indivual whose covariates are `covariates` and random effects `eta_i`.
#'
#' @details
#' The models used for the parameters are :
#' \deqn{ h_l(\psi_{li}) = h_l(\psi_{lpop})+X_i\beta_l + \eta_{li} }
#' with \eqn{h_l} the transformation, \eqn{\beta_l} the vector of covariates effect and  with \eqn{\eta_i} the random effects associated \eqn{\psi_l} parameter ;
#' \deqn{ g_k(\phi_{ki}) = g_k(\phi_{kpop})+X_i \gamma_l + \eta_{ki}}
#' with \eqn{g_k} the transformation and  \eqn{\gamma_k} the vector of covariates effect associated \eqn{\phi_k} parameter.
#'
#' @param theta model parameter contain phi_pop, psi_pop, gamma, beta ; list with name `phi_pop`, `psi_pop`, `gamma`, `beta` ;
#'
#' - `phi_pop` named vector of population parameters without r.e ;
#'
#' - `psi_pop` named vector of population parameters with r.e ;
#'
#' - `gamma` named list of vector of covariates effects for  `phi_pop` parameters, if NULL no covariates effect on parameters. ;
#'
#' - `beta` named list of vector of covariates effects for each `psi_pop`, if NULL no covariates effect on parameters.
#' @param covariates line matrix/data.frame of individual covariates ;
#' @param eta_i named vector of random effect for each `psi` parameter ;
#' @param transfo named list of transformation for `phi_pop` and `psi_pop` parameters ; if NULL id transformation is used;
#' @param transfo.inv named list of inverse transformation for `phi_pop` and `psi_pop` parameters ; if NULL id transformation is used.
#'
#' @return a list with `phi_i` and `psi_i` parameters.
#' @export
#' @seealso [Clairon()], [Pasin()]
#'
#' @examples
#'
#' phi_pop = c(delta_S = 0.231, delta_L = 0.000316)
#' psi_pop = c(delta_Ab = 0.025,phi_S = 3057, phi_L = 16.6)
#' gamma = NULL
#' covariates = data.frame(cAGE = runif(1,-15,15), G1 = rnorm(1), G2 = rnorm(1))
#' beta = list(delta_Ab=c(0,1.2,0),phi_S = c(0.93,0,0),phi_L=c(0,0,0.8))
#'
#' theta=list(phi_pop = phi_pop,psi_pop = psi_pop,gamma = gamma, beta = beta)
#' eta_i = c(delta_Ab = rnorm(1,0,0.3),phi_S=rnorm(1,0,0.92),phi_L=rnorm(1,0,0.85))
#' transfo = list(delta_Ab=log,phi_S=log,phi_L=log)
#' transfo.inv = list(delta_Ab = exp,phi_S=exp,phi_L=exp)
#'
#' indParm(theta,covariates,eta_i,transfo,transfo.inv)

indParm <- function(theta,covariates,eta_i,transfo,transfo.inv){
  # should add check parameter : a checker : theta names / beta vs covariates format / eta_i vs psi_pop format vs transfo vs transfo.inv format + names

  phi_i = numeric()
  psi_i = numeric()

  phi.names = names(theta$phi_pop)
  psi.names = names(theta$psi_pop)

  if(!is.matrix(covariates)){
    covariates <- as.matrix(covariates)
  }

  ## Compute phi parameters -> parameters with no r.e
  for(p in phi.names){
    trans = transfo[[p]]
    if(is.null(trans)){
      trans <- trans.inv <- function(x){x}
    }else{trans.inv <- transfo.inv[[p]]}
    if(is.null(theta$gamma[[p]])){
      phi_i[[p]] <- trans.inv(trans(theta$phi_pop[[p]]))
    }else{
      phi_i[[p]] <- trans.inv(trans(theta$phi_pop[[p]])+as.numeric(covariates%*%matrix(theta$gamma[[p]],ncol=1)))
    }
  }

  ## Compute psi parameters -> parameters with r.e.
  for(p in psi.names){
    trans = transfo[[p]]
    if(is.null(trans)){
      trans <- trans.inv <- function(x){x}
    }else{trans.inv <- transfo.inv[[p]]}
    if(is.null(theta$beta[[p]])){
      psi_i[[p]] <- trans.inv(trans(theta$psi_pop[[p]])+eta_i[[p]])
    }else{
      psi_i[[p]] <- trans.inv(trans(theta$psi_pop[[p]])+as.numeric(covariates%*%matrix(theta$beta[[p]],ncol=1)) + eta_i[[p]])
    }
  }
  return(list(phi_i=phi_i,psi_i=psi_i))
}

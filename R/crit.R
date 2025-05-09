#' BIC for remix object
#'
#' Computes bayesian information criterion from the output of \code{\link{remix}} as
#' \deqn{BIC = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+\log(N)P}
#' where \eqn{P} is the total number of parameters estimated, \eqn{N} the number of subject and \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model.
#'
#' @param object output of \code{\link{remix}}.
#' @param ... additional arguments.
#'
#' @references Schwarz, G. 1978. Estimating the dimension of a model. The annals of statistics 6 (2): 461–464
#' @returns BIC.
#' @export
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
#' res = remix(project = project,
#'             dynFUN = dynFUN_demo,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             eps2=1,
#'             lambda=lambda)
#'
#' BIC(res)
#' }
BIC.remix <- function(object,...){
  P = sum(object$finalRes$alpha!=0) + length(object$info$param.toprint)
  return(-2*object$finalRes$LL+log(object$info$N)*P)
}


#' AIC for remix object
#'
#' Computes akaike information criterion from the output of \code{\link{remix}} as
#' \deqn{AIC = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+k\times P}
#' where \eqn{P} is the total number of parameters estimated and \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model.
#'
#' @param object output of \code{\link{remix}}.
#' @param ... additional arguments.
#' @param k 	numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#'
#' @references Akaike, H. 1998. Information theory and an extension of the maximum likelihood principle, Selected papers of hirotugu akaike, 199–213. New York: Springer.
#' @returns AIC.
#' @export
#'
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
#' res = remix(project = project,
#'             dynFUN = dynFUN_demo,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             eps2=1,
#'             lambda=lambda)
#'
#' AIC(res)
#' }
AIC.remix <- function(object,...,k){
  P = sum(object$finalRes$alpha!=0) + length(object$info$param.toprint)
  return(-2*object$finalRes$LL+k*P)
}

#' eBIC
#'
#' Computes extended bayesian information criterion  as
#'  \deqn{ eBIC = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+P\log(N)+2\gamma\log(\binom(k,K))}
#' where \eqn{P} is the total number of parameters estimated, \eqn{N} the number of subject, \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model, \eqn{K} the number of submodel to explore (here the numbre of biomarkers tested) and \eqn{k} the numbre of biomarkers selected in the model.
#'
#' @param object output of \code{\link{remix}} or \code{\link{cv.remix}}.
#' @param ... opptional additional arguments.
#'
#' @references  Chen, J. and Z. Chen. 2008. Extended Bayesian information criteria for model selection with large model spaces. Biometrika 95 (3): 759–771
#' @returns eBIC.
#' @export
#'
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
#' res = remix(project = project,
#'             dynFUN = dynFUN_demo,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             eps2=1,
#'             lambda=lambda)
#'
#' eBIC(res)
#' }
eBIC <- function(object, ...) {
  UseMethod("eBIC", object)
}

#' @export
eBIC.remix <- function(object,gamma=1,...){
  P = length(object$info$param.toprint) + sum(object$finalRes$alpha!=0)
  return(-2*object$finalRes$LL+log(object$info$N)*P+2*gamma*log(choose(length(object$finalRes$alpha),sum(object$finalRes$alpha!=0))))
}

#' BICc
#'
#' Computes corrected bayesian information criterion  as
#'  \deqn{ BICc = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+P_R\log(N)+P_F\log(n_{tot})}
#' where \eqn{P_F} is the total number of parameters linked to fixed effects, \eqn{P_R} to random effects, \eqn{N} the number of subject, \eqn{n_tot} the total number of observations and \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model.
#'
#' @param object output of \code{\link{remix}} or \code{\link{cv.remix}}
#' @param ... opptional additional arguments.
#'
#' @references  Delattre M, Lavielle M, Poursat M-A. A note on BIC in mixed-effects models. Elect J Stat. 2014; 8(1): 456-475.
#' @returns BICc.
#' @export
#'
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
#' res = remix(project = project,
#'             dynFUN = dynFUN_demo,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             eps2=1,
#'             lambda=lambda)
#'
#' BICc(res)
#' }
BICc <- function(object, ...) {
  UseMethod("BICc", object)
}

#' @export
BICc.remix <- function(object,...){
  omega = object$info$param.toprint[stringr::str_detect(object$info$param.toprint,"omega_")]
  corr = object$info$param.toprint[stringr::str_detect(object$info$param.toprint,"corr_")]
  beta = object$info$param.toprint[stringr::str_detect(object$info$param.toprint,"beta_") & sub("^[^_]*_(.*)_[^_]*$", "\\1", object$info$param.toprint) %in% stringr::str_remove_all(omega,"omega_")]

  REP = Reduce(union,list(omega,corr,beta))
  FEP = setdiff(object$info$param.toprint,REP) # avec les alpha1 en moins car compté à part

  PF = length(FEP) + sum(object$finalRes$alpha!=0)
  PR = length(REP)
  return(-2*object$finalRes$LL+log(object$info$N)*PR+log(object$info$ntot)*PF)
}


# cv.REMIX-----------------------------------------------------------------

#' @export
BIC.cvRemix <- function(object,...){
  P = sapply(object$res,FUN=function(object){sum(object$alpha!=0)}) + length(object$info$param.toprint)
  return(-2*object$LL+log(object$info$N)*P)
}

#' @export
AIC.cvRemix <- function(object,...,k){
  P = sapply(object$res,FUN=function(object){sum(object$alpha!=0)}) + length(object$info$param.toprint)
  return(-2*object$LL+k*P)
}

#' @export
eBIC.cvRemix <- function(object,gamma=1,...){
  P =sapply(object$res,FUN=function(object){sum(object$alpha!=0)})+ length(object$info$param.toprint)
  return(-2*object$LL+log(object$info$N)*P+2*gamma*log(choose(length(object$info$alpha$alpha1),sapply(object$res,FUN=function(object){sum(object$alpha!=0)}))))
}

#' @export
BICc.cvRemix <- function(object,...){
  omega = object$info$param.toprint[stringr::str_detect(object$info$param.toprint,"omega_")]
  corr = object$info$param.toprint[stringr::str_detect(object$info$param.toprint,"corr_")]
  beta = object$info$param.toprint[stringr::str_detect(object$info$param.toprint,"beta_") & sub("^[^_]*_(.*)_[^_]*$", "\\1", object$info$param.toprint) %in% stringr::str_remove_all(omega,"omega_")]

  REP = Reduce(union,list(omega,corr,beta))
  FEP = setdiff(object$info$param.toprint,REP) # avec les alpha1 en moins car compté à part

  PF = length(FEP) + sapply(object$res,FUN=function(object){sum(object$alpha!=0)})
  PR = length(REP)

  return(-2*object$LL+log(object$info$N)*PR+log(object$info$ntot)*PF)
}

#'@export
logLik.remix <- function(object,...){
  return(object$finalRes$LL)
}

#' @export
logLik.cvRemix <- function(object,...){
  return(object$LL)
}

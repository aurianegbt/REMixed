#' BIC
#'
#' Computes bayesian information criterion from the output of \code{\link{remix}} as
#' \deqn{BIC = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+\log(N)P}
#' where \eqn{P} is the total number of parameters estimated, \eqn{N} the number of subject and \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model.
#'
#' @param x output of \code{\link{remix}} or \code{\link{cv.remix}}.
#' @param ... opptional additional arguments.
#'
#' @references   Schwarz, G. 1978. Estimating the dimension of a model. The annals of statistics 6 (2): 461–464
#' @returns BIC.
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
#' BIC(res)
#' }
BIC <- function(x, ...) {
  UseMethod("BIC", x)
}

#' @export
BIC.remix <- function(x,...){
  P = sum(x$finalRes$alpha!=0) + length(x$info$param.toprint)
  return(-2*x$finalRes$LL+log(x$info$N)*P)
}

#' @export
AIC.remix <- function(x,...){
  P = sum(x$finalRes$alpha!=0) + length(x$info$param.toprint)
  return(-2*x$finalRes$LL+2*P)
}

#' eBIC
#'
#' Computes extended bayesian information criterion  as
#'  \deqn{ eBIC = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+P\log(N)+2\gamma\log(\binom(k,K))}
#' where \eqn{P} is the total number of parameters estimated, \eqn{N} the number of subject, \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model, \eqn{K} the number of submodel to explore (here the numbre of biomarkers tested) and \eqn{k} the numbre of biomarkers selected in the model.
#'
#' @param x output of \code{\link{remix}} or \code{\link{cv.remix}}.
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
eBIC <- function(x, ...) {
  UseMethod("eBIC", x)
}

#' @export
eBIC.remix <- function(x,gamma=1,...){
  P = length(x$info$param.toprint) + sum(x$finalRes$alpha!=0)
  return(-2*x$finalRes$LL+log(x$info$N)*P+2*gamma*log(choose(length(x$finalRes$alpha),sum(x$finalRes$alpha!=0))))
}

#' BICc
#'
#' Computes corrected bayesian information criterion  as
#'  \deqn{ BICc = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+P_R\log(N)+P_F\log(n_{tot})}
#' where \eqn{P_F} is the total number of parameters linked to fixed effects, \eqn{P_R} to random effects, \eqn{N} the number of subject, \eqn{n_tot} the total number of observations and \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model.
#'
#' @param x output of \code{\link{remix}} or \code{\link{cv.remix}}
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
BICc <- function(x, ...) {
  UseMethod("BICc", x)
}

#' @export
BICc.remix <- function(x,...){
  omega = x$info$param.toprint[stringr::str_detect(x$info$param.toprint,"omega_")]
  corr = x$info$param.toprint[stringr::str_detect(x$info$param.toprint,"corr_")]
  beta = x$info$param.toprint[stringr::str_detect(x$info$param.toprint,"beta_") & sub("^[^_]*_(.*)_[^_]*$", "\\1", x$info$param.toprint) %in% stringr::str_remove_all(omega,"omega_")]

  REP = Reduce(union,list(omega,corr,beta))
  FEP = setdiff(x$info$param.toprint,REP) # avec les alpha1 en moins car compté à part

  PF = length(FEP) + sum(x$finalRes$alpha!=0)
  PR = length(REP)
  return(-2*x$finalRes$LL+log(x$info$N)*PR+log(x$info$ntot)*PF)
}


# cv.REMIX-----------------------------------------------------------------

#' @export
BIC.cvRemix <- function(x,...){
  P = sapply(x$res,FUN=function(x){sum(x$alpha!=0)}) + length(x$info$param.toprint)
  return(-2*x$LL+log(x$info$N)*P)
}

#' @export
AIC.cvRemix <- function(x,...){
  P = sapply(x$res,FUN=function(x){sum(x$alpha!=0)}) + length(x$info$param.toprint)
  return(-2*x$LL+2*P)
}

#' @export
eBIC.cvRemix <- function(x,gamma=1,...){
  P =sapply(x$res,FUN=function(x){sum(x$alpha!=0)})+ length(x$info$param.toprint)
  return(-2*x$LL+log(x$info$N)*P+2*gamma*log(choose(length(x$info$alpha$alpha1),sapply(x$res,FUN=function(x){sum(x$alpha!=0)}))))
}

#' @export
BICc.cvRemix <- function(x,...){
  omega = x$info$param.toprint[stringr::str_detect(x$info$param.toprint,"omega_")]
  corr = x$info$param.toprint[stringr::str_detect(x$info$param.toprint,"corr_")]
  beta = x$info$param.toprint[stringr::str_detect(x$info$param.toprint,"beta_") & sub("^[^_]*_(.*)_[^_]*$", "\\1", x$info$param.toprint) %in% stringr::str_remove_all(omega,"omega_")]

  REP = Reduce(union,list(omega,corr,beta))
  FEP = setdiff(x$info$param.toprint,REP) # avec les alpha1 en moins car compté à part

  PF = length(FEP) + sapply(x$res,FUN=function(x){sum(x$alpha!=0)})
  PR = length(REP)

  return(-2*x$LL+log(x$info$N)*PR+log(x$info$ntot)*PF)
}

#'@export
logLik.remix <- function(x,...){
  return(x$finalRes$LL)
}

#' @export
logLik.cvRemix <- function(x,...){
  return(x$LL)
}

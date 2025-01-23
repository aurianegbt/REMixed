#' BIC for remix object
#'
#' Computes bayesian information criterion from the output of \code{\link{remix}}.
#'
#' @param x output of \code{\link{remix}}.
#'
#' @returns BIC.
#' @export
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
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
BIC.remix <- function(x){
  return(-2*x$finalRes$LL+log(x$info$N)*sum(x$finalRes$alpha!=0))
}

#' AIC for remix object
#'
#' Computes akaike information criterion from the output of \code{\link{remix}}.
#'
#' @param x output of \code{\link{remix}}.
#'
#' @returns AIC.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
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
AIC.remix <- function(x){
  return(-2*x$finalRes$LL+2*sum(x$finalRes$alpha!=0))
}

#' eBIC
#'
#' Computes extended bayesian information criterion .
#'
#' @param x output of \code{\link{remix}} or \code{\link{cv.remix}}.
#'
#' @returns eBIC.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
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

#' eBIC for remix object
#'
#' Computes extended bayesian information criterion from the output of \code{\link{remix}}.
#'
#' @param x output of \code{\link{remix}}.
#'
#' @returns eBIC.
#' @references Chen, J., & Chen, Z. (2008). Extended Bayesian Information Criteria for Model Selection with Large Model Spaces. Biometrika, 95(3), 759–771. http://www.jstor.org/stable/20441500
#' @export
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
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
eBIC.remix <- function(x,gamma=1){
  return(-2*x$finalRes$LL+log(x$info$N)*sum(x$finalRes$alpha!=0)+2*gamma*log(choose(length(x$finalRes$alpha),sum(x$finalRes$alpha!=0))))
}

#' AICc
#'
#' Computes corrected akaike information criterion .
#'
#' @param x output of \code{\link{remix}} or \code{\link{cv.remix}}.
#'
#' @returns AICc.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
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
#' AICc(res)
#' }
AICc <- function(x, ...) {
  UseMethod("AICc", x)
}

#' AICc for remix object
#'
#' Computes corrected akaike information criterion from the output of \code{\link{remix}}.
#'
#' @param x output of \code{\link{remix}}.
#'
#' @returns AICc.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
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
#' AICc(res)
#' }
AICc.remix <- function(x){
  k=sum(x$finalRes$alpha!=0)
  return(-2*x$finalRes$LL+2*k+2*(k*(k+1))/(x$info$N-k-1))
}


# cv.REMIX-----------------------------------------------------------------

#' BIC for cvRemix object
#'
#' Computes bayesian information criterion from the output of \code{\link{cv.rremix}}.
#'
#' @param x output of \code{\link{remix}}.
#'
#' @returns BIC across the grid.
#' @export
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
#'
#' res = cv.remix(project = project,
#'                dynFUN = dynFUN_demo,
#'                y = y,
#'                ObsModel.transfo = ObsModel.transfo,
#'                alpha = alpha,
#'                selfInit = TRUE,
#'                eps1=10**(-2),
#'                ncores=8,
#'                nlambda=8,
#'                eps2=1)
#'
#' BIC(res)
#' }
BIC.cvRemix <- function(x){
  return(-2*x$LL+log(x$info$N)*sapply(x$res,FUN=function(x){sum(x$alpha!=0)}))
}

#' AIC for cvRemix object
#'
#' Computes akaike information criterion from the output of \code{\link{cv.remix}}.
#'
#' @param x output of \code{\link{remix}}.
#'
#' @returns AIC across the grid.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
#'
#' res = cv.remix(project = project,
#'                dynFUN = dynFUN_demo,
#'                y = y,
#'                ObsModel.transfo = ObsModel.transfo,
#'                alpha = alpha,
#'                selfInit = TRUE,
#'                eps1=10**(-2),
#'                ncores=8,
#'                nlambda=8,
#'                eps2=1)
#'
#' AIC(res)
#' }
AIC.cvRemix <- function(x){
  return(-2*x$LL+2*sapply(x$res,FUN=function(x){sum(x$alpha!=0)}))
}

#' eBIC for cvRemix object
#'
#' Computes extended bayesian information criterion from the output of \code{\link{cv.remix}}.
#'
#' @param x output of \code{\link{remix}}.
#'
#' @returns eBIC accross the grid.
#' @references Chen, J., & Chen, Z. (2008). Extended Bayesian Information Criteria for Model Selection with Large Model Spaces. Biometrika, 95(3), 759–771. http://www.jstor.org/stable/20441500
#' @export
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
#'
#' res = cv.remix(project = project,
#'                dynFUN = dynFUN_demo,
#'                y = y,
#'                ObsModel.transfo = ObsModel.transfo,
#'                alpha = alpha,
#'                selfInit = TRUE,
#'                eps1=10**(-2),
#'                ncores=8,
#'                nlambda=8,
#'                eps2=1)
#'
#' eBIC(res)
#' }
eBIC.cvRemix <- function(x,gamma=1){
  return(-2*x$LL+log(x$info$N)*sapply(x$res,FUN=function(x){sum(x$alpha!=0)})+2*gamma*log(choose(length(x$info$alpha$alpha1),sapply(x$res,FUN=function(x){sum(x$alpha!=0)}))))
}

#' AICc for cvRemix object
#'
#' Computes corrected akaike information criterion from the output of \code{\link{cv.remix}}.
#'
#' @param x output of \code{\link{cv.remix}}.
#'
#' @returns AICc accross the grid.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(Ab=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,Ab=1000)
#'
#' res = cv.remix(project = project,
#'                dynFUN = dynFUN_demo,
#'                y = y,
#'                ObsModel.transfo = ObsModel.transfo,
#'                alpha = alpha,
#'                selfInit = TRUE,
#'                eps1=10**(-2),
#'                ncores=8,
#'                nlambda=8,
#'                eps2=1)
#'
#' AICc(res)
#' }
AICc.cvRemix <- function(x){
  k=sapply(x$res,FUN=function(x){sum(x$alpha!=0)})
  return(-2*x$LL+2*k+2*(k*(k+1))/(x$info$N-k-1))
}


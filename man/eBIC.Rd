% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crit.R
\name{eBIC}
\alias{eBIC}
\title{eBIC}
\usage{
eBIC(object, ...)
}
\arguments{
\item{object}{output of \code{\link{remix}} or \code{\link{cv.remix}}.}

\item{...}{opptional additional arguments.}
}
\value{
eBIC.
}
\description{
Computes extended bayesian information criterion  as
 \deqn{ eBIC = -2\mathcal{LL}_{y}(\hat\theta,\hat\alpha)+P\log(N)+2\gamma\log(\binom(k,K))}
where \eqn{P} is the total number of parameters estimated, \eqn{N} the number of subject, \eqn{\mathcal{LL}_{y}(\hat\theta,\hat\alpha)} the log-likelihood of the model, \eqn{K} the number of submodel to explore (here the numbre of biomarkers tested) and \eqn{k} the numbre of biomarkers selected in the model.
}
\examples{
\dontrun{
project <- getMLXdir()

ObsModel.transfo = list(S=list(AB=log10),
                        linkS="yAB",
                        R=rep(list(S=function(x){x}),5),
                        linkR = paste0("yG",1:5))

alpha=list(alpha0=NULL,
           alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))

y = c(S=5,AB=1000)
lambda = 1440

res = remix(project = project,
            dynFUN = dynFUN_demo,
            y = y,
            ObsModel.transfo = ObsModel.transfo,
            alpha = alpha,
            selfInit = TRUE,
            eps1=10**(-2),
            eps2=1,
            lambda=lambda)

eBIC(res)
}
}
\references{
Chen, J. and Z. Chen. 2008. Extended Bayesian information criteria for model selection with large model spaces. Biometrika 95 (3): 759-771.
}

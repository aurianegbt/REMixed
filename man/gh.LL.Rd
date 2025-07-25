% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gh_LL.R
\name{gh.LL}
\alias{gh.LL}
\title{Adaptive Gauss-Hermite approximation of log-likelihood derivatives}
\usage{
gh.LL(
  dynFUN,
  y,
  mu = NULL,
  Omega = NULL,
  theta = NULL,
  alpha1 = NULL,
  covariates = NULL,
  ParModel.transfo = NULL,
  ParModel.transfo.inv = NULL,
  Sobs = NULL,
  Robs = NULL,
  Serr = NULL,
  Rerr = NULL,
  ObsModel.transfo = NULL,
  data = NULL,
  n = NULL,
  prune = NULL,
  parallel = TRUE,
  ncores = NULL,
  onlyLL = FALSE,
  verbose = TRUE,
  precBits = 10
)
}
\arguments{
\item{dynFUN}{function computing the dynamics of interest for a set of parameters. This function need to contain every sub-function that it may needs (as it is called in a foreach loop). The output of this function need to return a data.frame with \code{time} : as first columns and named dynamics in other columns. It must take in input : \itemize{\item\code{y} : a named vector with the initial condition. The names are the dynamics names.
\item\code{parms} : a named vector of parameter.
\item\code{time} : vector a timepoint.}

See \code{\link{dynFUN_demo}}, \code{\link{model.clairon}}, \code{\link{model.pasin}} or \code{\link{model.pk}} for examples.}

\item{y}{initial condition of the mechanism model, conform to what is asked in \code{dynFUN}. If regressor used in Monolix provided a named list of vector of individual initial conditions. Each vector need to be of length 1 (same for all), or exactly the numbre of individuals (range in the same order as their id).}

\item{mu}{list of individuals random effects estimation (vector of r.e. need to be named by the parameter names), use to locate the density mass; (optional, see description).}

\item{Omega}{list of individuals estimated standard deviation diagonal matrix (matrix need to have rows and columns named by the parameter names), use to locate the density mass; (optional, see description).}

\item{theta}{list of model parameters containing (see details) \itemize{\item\code{phi_pop} : named vector with the population parameters with no r.e. \eqn{(\phi_{l\ pop})_{l\leq L}} (NULL if none) ;
\item\code{psi_pop} : named vector with the population parameters with r.e. \eqn{(\psi_{l\ pop})_{l\leq m}} ;
\item\code{gamma} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with no r.e. ;
\item\code{beta} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with r.e..
\item\code{alpha0} : named vector of \eqn{(\alpha_{0k})_{k\leq K}} parameters (names are identifier of the observation model, such as in a Monolix project);
\item\code{omega} : named vector of estimated r.e. standard deviation;}
(optional, see description).}

\item{alpha1}{named vector of regulatization parameters \eqn{(\alpha_{1k})_{k\leq K}}, with identifier of observation model as names, (optional, see description).}

\item{covariates}{matrix of individual covariates (size N x n). Individuals must be sorted in the same order than in \code{mu} and \code{Omega}, (optional, see description).}

\item{ParModel.transfo}{named list of transformation functions \eqn{(h_l)_{l\leq m}} and \eqn{(s_k)_{k\leq K}} for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop}, missing entries are set by default to the identity function ; optional, see description).}

\item{ParModel.transfo.inv}{Named list of inverse transformation functions for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop} ; optional, see description).}

\item{Sobs}{list of individuals trajectories for the direct observation models \eqn{(Y_{pi})_{p \leq P,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{p\leq P} data.frame with time  \eqn{(t_{pij})_{j\leq n_{ip}}} and observations  \eqn{(Y_{pij})_{j\leq n_{ip}}}. Each data.frame is named with the observation model identifiers.}

\item{Robs}{list of individuals trajectories for the latent observation models \eqn{(Z_{ki})_{k \leq K,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{k\leq K} data.frame with time  \eqn{(t_{kij})_{j\leq n_{ik}}} and observations  \eqn{(Z_{kij})_{j\leq n_{ik}}}. Each data.frame is named with the observation model identifiers.}

\item{Serr}{named vector of the estimated error mocel constants \eqn{(\varsigma_p)_{p\leq P}} with observation model identifiers as names.}

\item{Rerr}{named vector of the estimated error mocel constants \eqn{(\sigma_k)_{k\leq K}} with observation model identifiers as names.}

\item{ObsModel.transfo}{list containing two lists of transformations and two vectors linking each transformations to their observation model name in the Monolix project. The list should include identity transformations and be named \code{S} and \code{R}. The two vectors should be named \code{linkS} and \code{linkR}.

Both \code{S} (for the direct observation models) and \code{linkS}, as well as \code{R} (for latent process models) and \code{linkR}, must have the same length.

\itemize{
  \item\code{S}: a list of transformations for the direct observation models. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic is observed (from \code{dynFUN});  \item\code{linkS} : a vector specifying the observation model names (that is used in the monolix project, \code{alpha1}, etc.) for each transformation, in the same order as in \code{S};

  \item\code{R}: similarly, a list of transformations for the latent process models. Although currently there is only one latent dynamic, each \eqn{s_k, k\leq K} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from \code{dynFUN}; \item \code{linkR} : a vector specifying the observation model names for each transformation, in the same order as in \code{R}.
}}

\item{data}{output from \code{\link{readMLX}} containing parameters "\code{mu}", "\code{Omega}", "\code{theta}", "\code{alpha1}", "\code{covariates}", "\code{ParModel.transfo}", "\code{ParModel.transfo.inv}", "\code{Sobs}", "\code{Robs}", "\code{Serr}", "\code{Rerr}", "\code{ObsModel.transfo}" extract from a monolix project.}

\item{n}{number of points per dimension to use for the Gauss-Hermite quadrature rule.}

\item{prune}{integer between 0 and 1, percentage of pruning for the Gauss-Hermite quadrature rule (default NULL).}

\item{parallel}{logical, if computation should be done in parallel.}

\item{ncores}{number of cores to use for parallelization, default will detect the number of cores available.}

\item{onlyLL}{logical, if only the log-likelihood should be computed (and not \eqn{\partial_{\alpha_1} LL} or \eqn{\partial_{\alpha_1}^2 LL}).}

\item{verbose}{logical, if progress bar should be printed through the computation.}

\item{precBits}{precision if needed}
}
\value{
A list with the approximation by Gauss-Hermite quadrature of the likelihood \code{L}, the log-likelihood \code{LL}, the gradient of the log-likelihood \code{dLL}, and the Hessian of the log-likelihood \code{ddLL} at the point \eqn{\theta, \alpha} provided.
}
\description{
Computes Adaptive Gauss-Hermite approximation of the log-likelihood and its derivatives in NLMEM with latent observation processes, see \code{\link{REMixed-package}} for details on the model.
}
\details{
Based on notation introduced \code{\link{REMixed-package}}. The log-likelihood of the model \eqn{LL(\theta,\alpha_1)} for a set of population parameters \eqn{\theta} and regulatization parameters \eqn{\alpha_1} is estimated using Adaptative Gausse-Hermite quadrature, using  conditional distribution estimation to locate the mass of the integrand. If the project has been initialized as a Monolix project, the user can use \code{\link{readMLX}} function to retrieve all the project information needed here.
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

data <- readMLX(project,ObsModel.transfo,alpha)

LL <- gh.LL(dynFUN = dynFUN_demo,
            y = c(S=5,AB=1000),
            ObsModel.transfo=ObsModel.transfo,
            data = data)

print(LL)
}
}

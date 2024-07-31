#' Adaptive Gauss-Hermite Approximation of Log-Likelihood Derivatives
#'
#' @description
#' Computes the adaptive Gauss-Hermite approximation of the log-likelihood and its derivatives in NLMEM with latent observation processes.
#'
#' @details
#' Suppose we have a differential system of equations containing variables \eqn{(S_{p})_{p\leq P}} and \eqn{R} that depend on some parameters. These dynamics are described by `dynFUN`. We model the process over time for \eqn{i \leq N} individuals, resulting from the differential system for a set of parameters \eqn{\phi_i} and \eqn{(\psi_{li})_{l\leq m, i\leq N}} for individual \eqn{i \leq N}, as \eqn{S_{p}(\cdot, \phi_i, (\psi_{li})_{l\leq m}) = S_{pi}(\cdot)}, \eqn{p \leq P}, and \eqn{R(\cdot, \phi_i, (\psi_{li})_{l\leq m}) = R_i(\cdot)}.
#'
#' Parameters are described as \deqn{h_l(\psi_{li}) = h_l(\psi_{lpop}) + X_i \beta_l + \eta_{li}}, where the covariates of individual \eqn{i \leq N} are contained in the rows of the matrix `covariates`, and random effects \deqn{\eta_i = (\eta_{li})_{l \leq m} \overset{iid}{\sim} \mathcal{N}(\mu_i, \Omega_i)} for \eqn{i \leq N}, where \eqn{\mu_i} (in the list `mu`) is the estimated random effect of individual \eqn{i} and \eqn{\Omega_i} (in the list `Omega`) is the diagonal matrix of estimated standard deviations of random effects for individual \eqn{i}. The population parameters \eqn{\psi_{pop} = (\psi_{lpop})_{l \leq m}} (contained in `theta` through `psi_pop`) and \eqn{\beta = (\beta_l)_{l \leq m}} (contained in `theta` through `beta`) form the vector of covariate effects on parameters.
#'
#' The remaining population parameters of the structural model, which do not have random effects, are denoted by \eqn{(\phi_i)_{i \leq N}}, and are defined as \eqn{\phi_i = \phi_{pop} + X_i \gamma}. These can depend on covariate effects or be constant across the population. The population parameters \eqn{\phi_{pop}} are contained in `theta` through `phi_pop`, and the covariate effects coefficients \eqn{\gamma} are contained in `theta` through `gamma`.
#'
#' We assume that individual trajectories \eqn{(S_{pi})_{p \leq P, i \leq N}} are observed through a direct observation model, up to a transformation \eqn{g_p}, \eqn{p \leq P}, at different times \eqn{(t_{pij})_{i \leq N, p \leq P, j \leq n_{ip}}}: \deqn{Y_{pij} = g_p(S_{pi}(t_{pij})) + \epsilon_{pij}}, contained in the list `Sobs` with error \eqn{\epsilon_p = (\epsilon_{pij}) \overset{iid}{\sim} \mathcal{N}(0, \varsigma_p^2)} for \eqn{p \leq P} (with \eqn{(\varsigma_p)_{p \leq P}} contained in `Serr`).
#'
#' The individual trajectory \eqn{(R_{i})_{i \leq N}} is observed through latent processes, up to a transformation \eqn{s_k}, \eqn{k \leq K}, observed at \eqn{(t_{kij})_{i \leq N, k \leq K, j \leq n_{kij}}}: \deqn{Z_{kij} = \alpha_{k0} + \alpha_{k1} s_k(R_i(t_{kij})) + \varepsilon_{kij}} (contained in the list `Robs`), where \eqn{\varepsilon_k \overset{iid}{\sim} \mathcal{N}(0, \sigma_k^2)} (with \eqn{(\sigma_k)_{k \leq K}} contained in `Rerr`).
#'
#' @param dynFUN Dynamic function. This function need to contain every sub-function that it may needs as it is called in a foreach loop. The output of this function need to return a data.frame with `time` as first columns and named dynamics in other columns. It must take in input : \itemize{\item `y` a named vector with the initial condition. The names are the dynamics names.
#' \item `parms` a named vector of parameter.
#' \item `time` vector a timepoint.}
#'
#' See \code{\link{dynFUN_demo}}, \code{\link{Clairon}}, \code{\link{Pasin}} or \code{\link{PK}} for example.
#' @param y Initial condition of the model, conforming to what is required in `dynFUN`.
#' @param theta Model parameters : contains `phi_pop` (size L), `psi_pop` (size m), `gamma` (list of size L, each containing a vector of size m), `beta` (list of size m, each containing a vector of size m), and `alpha0` (size K).
#' @param alpha1 Regularization parameter (size K). It should be named with observation model names.
#' @param ParModel.transfo Named list of transformation functions for the individual parameter model (names must be consistent with `phi_pop` and `psi_pop`, missing entries default to identity).
#' @param ParModel.transfo.inv Named list of inverse transformation functions for the individual parameter model (names must be consistent with `phi_pop` and `psi_pop`, missing entries default to identity).
#' @param Serr Named ector of size P containing estimated error model constants (must be in the same order as in `Sobs`). It should be named with observation model names.
#' @param Rerr Vector of size K containing estimated error model constants (must be in the same order as in `Robs`). It should be named with unique observation model names.
#' @param ObsModel.transfo List containing two lists of transformations and two vectors with their corresponding links to the observation models names. The list should include identity transformations and be named `S` and `R`. It should also include two vectors, `linkS` and `linkR`, which link to the observation models.
#'
#' Both `ObsModel.transfo$S` (for direct observation models) and `ObsModel.transfo$linkS`, as well as `ObsModel.transfo$R` (for latent process models) and `ObsModel.transfo$linkR`, must have the same length.
#'
#' \itemize{
#'   \item `ObsModel.transfo$S`: A list of transformations for the direct observation model. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic is observed. The `linkS` vector specifies the observation model name given in other inputs for each transformation, in the same order as in `ObsModel.transfo$S`.
#'
#'   \item `ObsModel.transfo$R`: A list of transformations for the latent process model. Although currently there is only one latent dynamic, each \eqn{s_k} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from `dynFUN`. The `linkR` vector specifies the observation model names for each transformation, in the same order as in `ObsModel.transfo$R`.
#' }
#' @param n Number of points for the adaptative Gaussian quadrature (default is \code{floor(100^(1/length(theta$psi_pop)))}).
#' @param prune Percentage for pruning in [0, 1] (default is NULL).
#' @param mu List of named vectors of individual random effects estimations (size N, each vector has size m). Individuals must be sorted in the same order than `Omega` and `covariates`.
#' @param Omega List of named diagonal matrices of individual standard deviation estimations (size N, each element is size m x m). Individuals must be sorted in the same order than `mu` and `covariates`.
#' @param covariates Matrix of individual covariates (size N x n). Individuals must be sorted in the same order than `mu` and `Omega`.
#' @param Sobs List (size N) of lists of direct observations (size P), each element containing time and observation (in column 3). The element list name must be the observation variable names.
#' @param Robs List (size N) of latent observations (size K), each element containing time and observation (in column 3). The element list name must be the observation variable names. (same as in `ObsModel.transfo$link[...]`, etc.)
#' @param ncores Number of cores for parallelization (default is NULL).
#' @param data A list containing all `mu`, `Omega`, `theta`, `alpha1`, `covariates`, `ParModel.transfo`, `ParModel.transfo.inv`, `Sobs`, `Robs`, `Serr`, `Rerr`, and `ObsModel.transfo`, which can be obtained using \code{\link{readMLX}} from a monolix project. If any of these arguments are specified, they will be used in priority.
#' @param onlyLL If only the log-likelihood needs to be returned (default is FALSE).
#' @param parallel If the computation should be done in parallel (default is TRUE).
#' @param verbose if TRUE, progress bar appears.
#'
#' @return A list with the approximation by Gauss-Hermite quadrature of the likelihood `L`, the log-likelihood `LL`, the gradient of the log-likelihood `dLL`, and the Hessian of the log-likelihood `ddLL` at the point \eqn{\theta, \alpha}.
#'
#' @export
#'
#' @examples
#' project <- getMLXdir()
#'
#'
#' ObsModel.transfo = list(S=list(AB=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' data <- readMLX(project,ObsModel.transfo,alpha)
#'
#' LL <- gh.LL(dynFUN = dynFUN_demo,
#'             y = c(S=5,AB=1000),
#'             ObsModel.transfo=ObsModel.transfo,
#'             data = data)
#'
#' print(LL)
gh.LL <- function(
    dynFUN,
    y,
    mu=NULL,
    Omega=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs=NULL,
    Robs=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    n=NULL,
    prune=NULL,
    parallel=TRUE,
    ncores=NULL,
    onlyLL=FALSE,
    verbose=TRUE){

  if(is.null(data)){
    test <- sapply(c("mu","Omega","theta","alpha1","covariates","ParModel.transfo","ParModel.transfo.inv","Sobs","Robs","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu) || is.null(Omega)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    for(d in 1:length(data)){
      test <- eval(parse(text=(paste0("is.null(",names(data)[d],")"))))
      if(test){
        eval(parse(text=paste0(names(data[d]),"<- data[[d]]")))
      }
    }
  }
  if(is.null(n)){
    n <- floor(100**(1/length(theta$psi_pop)))
  }

  if(parallel){
    if(!is.null(ncores)){
      cluster <- snow::makeCluster(ncores)
    }else{
      cluster <- snow::makeCluster(parallel::detectCores())
    }
    doSNOW::registerDoSNOW(cluster)
  }

  N=length(mu)

  i = 1

  ntasks <- N
  if(verbose){
    pb <- utils::txtProgressBar(max = ntasks, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{
    opts <- NULL
  }


  res = foreach::foreach(i = 1:N,.packages = "REMix",.options.snow=opts)%dopar%{
    if(0 %in% diag(Omega[[i]])){
      diag(Omega[[i]])[diag(Omega[[i]])==0] <- 10**(-5)
    }
    LLi = gh.LL.ind(mu_i = mu[[i]],
                    Omega_i = Omega[[i]],
                    theta = theta,
                    alpha1 = alpha1,
                    dynFUN = dynFUN,
                    y = y,
                    covariates_i = covariates[i,,drop=F],
                    ParModel.transfo = ParModel.transfo,
                    ParModel.transfo.inv = ParModel.transfo.inv,
                    Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==i,]}),
                    Robs_i = lapply(Robs,FUN=function(R){R[R$id==i,]}),
                    Serr = Serr,
                    Rerr = Rerr,
                    ObsModel.transfo = ObsModel.transfo,
                    n = n,
                    prune = prune,
                    onlyLL = onlyLL)


    return(LLi)
  }

  L = prod(sapply(res,FUN=function(ri){ri$Li}))

  LL = sum(sapply(res,FUN=function(ri){ri$LLi}))


  if(verbose)
    close(pb)
  if(parallel){
    snow::stopCluster(cluster)
  }

  if(!onlyLL){
    if(length(alpha1)!=1){
      dLL = rowSums(sapply(res,FUN=function(ri){ri$dLLi}))
      ddLL = matrix(rowSums(sapply(res,FUN=function(ri){ri$ddLLi})),ncol=length(dLL))
    }else{
      dLL = sum(sapply(res,FUN=function(ri){ri$dLLi}))
      ddLL = matrix(sum(sapply(res,FUN=function(ri){ri$ddLLi})),ncol=length(dLL))
    }
    return(list(L=L,LL=LL,dLL=dLL,ddLL=ddLL))
  }else{
    return(LL)
  }
}

gh.LL.ind <- function(
    dynFUN,
    y,
    mu_i=NULL,
    Omega_i=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates_i=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs_i=NULL,
    Robs_i=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    ind = NULL,
    n = NULL,
    prune=NULL,
    onlyLL=FALSE){

  if(is.null(data)){
    test <- sapply(c("mu_i","Omega_i","theta","alpha1","covariates_i","ParModel.transfo","ParModel.transfo.inv","Sobs_i","Robs_i","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu_i) || is.null(Omega_i)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    argsNONind = c("theta","alpha1","ParModel.transfo","ParModel.transfo.inv","Serr","Rerr","ObsModel.transfo")
    for(d in 1:length(argsNONind)){
      test <- eval(parse(text=paste0("is.null(",argsNONind[d],")")))
      if(test){
        eval(parse(text=paste0(argsNONind[d],"<- data[[argsNONind[d]]]")))
      }
    }
    argsInd = c("mu","Omega","Robs","Sobs","covariates")
    for(d in 1:length(argsInd)){
      test <- eval(parse(text=paste0("is.null(",paste0(argsInd[d],"_i"),")")))
      if(test){
        eval(parse(text=paste0(argsInd[d],"<- data[[argsInd[d]]]")))

      }
    }
    if((is.null(mu_i) || is.null(Omega_i) || is.null(Robs_i) || is.null(Sobs_i) || is.null(covariates_i)) & is.null(ind)){
      stop("Please provide individual information (mu_i,Omega_i,Sobs_i,Robs_i,covariates_i), or the individual id ind.")
    }
    if(is.null(mu_i)){
      mu_i = mu[[ind]]
    }
    if(is.null(Omega_i)){
      Omega_i = Omega[[ind]]
    }
    if(is.null(Sobs_i)){
      Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==ind,]})
    }
    if(is.null(Robs_i)){
      Robs_i = lapply(Robs,FUN=function(R){R[R$id==ind,]})
    }
    if(is.null(covariates_i)){
      covariates_i = covariates[ind,,drop=F]
    }
  }

  if(is.null(n)){
    n <- floor(100**(1/length(theta$psi_pop)))
  }

  dm = length(mu_i)

  if(dm!=1){
    mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
    detsq.omega = prod(theta$omega)
    root.omega = diag(1/theta$omega**2)
  }else{
    mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
    detsq.omega = as.numeric(theta$omega)
    root.omega = 1/as.numeric(theta$omega)**2
  }

  nd = length(mh.parm$Weights)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

  # Need to compute, for each eta_i x individual i, the dynamics of the model
  # eta_i = split(mh.parm$Points,1:nd)[[1]]
  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
    dyn_eta_i <- dynFUN(all.tobs,y,unlist(unname(PSI_i)))

    return(dyn_eta_i)
  }),paste0("eta_",1:nd))

  # Compute marginal latent density for each individual and value of RE
  R.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
    dyn.ei = dyn[[paste0("eta_",ei)]]


    res = prod(sapply(1:R.sz,FUN=function(k){
      yGk = ObsModel.transfo$linkR[k]
      sig = Rerr[[yGk]]
      tki = Robs_i[[yGk]]$time
      Zki = Robs_i[[yGk]][,yGk]

      a0 = theta$alpha0[[yGk]]
      a1 = alpha1[[yGk]]
      trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
      Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

      prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2))
    }))
  }),paste0("eta_",1:nd))

  # Compute marginal density for each individual and value of RE
  if(S.sz!=0){
    S.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
      eta_i = split(mh.parm$Points,1:nd)[[ei]]
      dyn.ei = dyn[[paste0("eta_",ei)]]

      res = prod(sapply(1:S.sz,FUN=function(p){
        Yp = ObsModel.transfo$linkS[p]
        sig = Serr[[Yp]]
        tpi = Sobs_i[[Yp]]$time
        Ypi = Sobs_i[[Yp]][,Yp]

        trs = ObsModel.transfo$S[which(ObsModel.transfo$linkS==Yp)]
        Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])

        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
      }))*(1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }else{
    S.margDensity = setNames(sapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
      (1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }

  # Compute individual log-Likelihood
  Li = sum(mh.parm$Weights*R.margDensity*S.margDensity)
  LLi = log(Li)

  if(!onlyLL){
    # Compute  gradient of individual log-Likelihood
    # dim.alpha = length(alpha1) # = K = length(Robs_i )
    dfact = lapply(1:R.sz,FUN=function(k){
      setNames(sapply(1:nd,FUN=function(ei){
        dyn.ei = dyn[[paste0("eta_",ei)]]
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        a1 = alpha1[[yGk]]
        trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

        return(1/(sig**2)*(sum(Rki*(Zki-a0-a1*Rki))))
      }),paste0("eta_",1:nd))
    })

    dLLi = sapply(dfact,FUN=function(f){
      return(1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
    })

    # Compute hessienne of individual log-Likelihood
    ddLLi = matrix(0,ncol=R.sz,nrow=R.sz)
    # diagonal term
    ddfact = lapply(1:R.sz,FUN=function(k){
      setNames(sapply(1:nd,FUN=function(ei){
        dyn.ei = dyn[[paste0("eta_",ei)]]
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        a1 = alpha1[[yGk]]
        trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

        return(
          (1/(sig**2)*(sum(Rki*(Zki-a0-a1*Rki))))**2 - 1/(sig**2)*(sum(Rki**2))
        )
      }),paste0("eta_",1:nd))
    })

    diag(ddLLi) = sapply(ddfact,FUN=function(f){
      return(1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
    })-(dLLi)**2

    # non diagonal term
    for(k1 in 1:R.sz){
      for(k2 in 1:R.sz){
        if(k1<k2){
          ddfact = setNames(sapply(1:nd,FUN=function(ei){
            dyn.ei = dyn[[paste0("eta_",ei)]]
            prod(sapply(c(k1,k2),FUN=function(k){
              yGk = ObsModel.transfo$linkR[k]

              sig = Rerr[[yGk]]
              tki = Robs_i[[yGk]]$time
              Zki = Robs_i[[yGk]][,yGk]

              a0 = theta$alpha0[[yGk]]
              a1 = alpha1[[yGk]]
              trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
              Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

              return(
                1/sig**2*(sum(Rki*(Zki-a0-a1*Rki)))
              )
            }))
          }),paste0("eta_",1:nd))

          ddLLi[k1,k2] <- ddLLi[k2,k1] <- 1/Li*sum(mh.parm$Weights*ddfact*R.margDensity*S.margDensity)-prod(dLLi[c(k1,k2)])
        }
      }
    }
  }else{
    dLLi <- ddLLi <- NULL
  }
  return(list(Li=Li,LLi=LLi,dLLi=dLLi,ddLLi=ddLLi))
}


# One dimensional Adaptative Hermite Gauss  ------------------------------------------
agauss.hermite <- function(n,mu=0,sd=1){
  gh <- fastGHQuad::gaussHermiteData(n)
  names(gh) <- c("Points","Weights")

  ww <- gh$Weights * sqrt(2) * sd * exp(gh$Points**2)
  xx <- mu + sd * sqrt(2) * gh$Points

  return(data.frame(Points=xx,Weights=ww))
}

# Multi-dimensional Adaptative Hermite Gauss -----------------------------------------
amgauss.hermite <- function(n,mu=rep(0,ncol(Omega)),Omega=diag(rep(1,length(mu))),prune=NULL){
  dm = length(mu)

  gh  <- as.data.frame(fastGHQuad::gaussHermiteData(n))

  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  # plot(pts, cex=-5/log(wts), pch=19)


  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }

  eig <- eigen(Omega)
  S = eig$vectors
  L = diag(sqrt(eig$values))
  rot <- S %*% L

  ppts  <-t( sqrt(2) * rot %*% t(pts) + mu)
  wwts  <- wts * sqrt(2)**dm * abs(det(rot)) * exp(apply(pts,1,FUN=function(x){t(x)%*%x}))


  return(list(Points=ppts, Weights=wwts))
}

# individual contribution to log-lik gradient in alpha1=0 ----------------------------
lambda.max.ind <- function(
    dynFUN,
    y,
    mu_i=NULL,
    Omega_i=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates_i=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs_i=NULL,
    Robs_i=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    ind = NULL,
    n = NULL,
    prune=NULL){

  if(is.null(data)){
    test <- sapply(c("mu_i","Omega_i","theta","alpha1","covariates_i","ParModel.transfo","ParModel.transfo.inv","Sobs_i","Robs_i","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu_i) || is.null(Omega_i)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    argsNONind = c("theta","alpha1","ParModel.transfo","ParModel.transfo.inv","Serr","Rerr","ObsModel.transfo")
    for(d in 1:length(argsNONind)){
      test <- eval(parse(text=paste0("is.null(",argsNONind[d],")")))
      if(test){
        eval(parse(text=paste0(argsNONind[d],"<- data[[argsNONind[d]]]")))
      }
    }
    argsInd = c("mu","Omega","Robs","Sobs","covariates")
    for(d in 1:length(argsInd)){
      test <- eval(parse(text=paste0("is.null(",paste0(argsInd[d],"_i"),")")))
      if(test){
        eval(parse(text=paste0(argsInd[d],"<- data[[argsInd[d]]]")))

      }
    }
    if((is.null(mu_i) || is.null(Omega_i) || is.null(Robs_i) || is.null(Sobs_i) || is.null(covariates_i)) & is.null(ind)){
      stop("Please provide individual information (mu_i,Omega_i,Sobs_i,Robs_i,covariates_i), or the individual id ind.")
    }
    if(is.null(mu_i)){
      mu_i = mu[[ind]]
    }
    if(is.null(Omega_i)){
      Omega_i = Omega[[ind]]
    }
    if(is.null(Sobs_i)){
      Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==ind,]})
    }
    if(is.null(Robs_i)){
      Robs_i = lapply(Robs,FUN=function(R){R[R$id==ind,]})
    }
    if(is.null(covariates_i)){
      covariates_i = covariates[ind,,drop=F]
    }
  }

  if(is.null(n)){
    n <- floor(100**(1/length(theta$psi_pop)))
  }

  dm = length(mu_i)

  if(dm!=1){
    mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
    detsq.omega = prod(theta$omega)
    root.omega = diag(1/theta$omega**2)
  }else{
    mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
    detsq.omega = as.numeric(theta$omega)
    root.omega = 1/as.numeric(theta$omega)**2
  }

  nd = length(mh.parm$Weights)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
    dyn_eta_i <- dynFUN(all.tobs,y,unlist(unname(PSI_i)))

    return(dyn_eta_i)
  }),paste0("eta_",1:nd))

  # Compute marginal latent density for each individual and value of RE
  R.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
    dyn.ei = dyn[[paste0("eta_",ei)]]


    res = prod(sapply(1:R.sz,FUN=function(k){
      yGk = ObsModel.transfo$linkR[k]
      sig = Rerr[[yGk]]
      tki = Robs_i[[yGk]]$time
      Zki = Robs_i[[yGk]][,yGk]

      a0 = theta$alpha0[[yGk]]
      a1 = alpha1[[yGk]]
      trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
      Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

      prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2))
    }))
  }),paste0("eta_",1:nd))

  # A(Y|theta0)
  if(S.sz!=0){
    S.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
      eta_i = split(mh.parm$Points,1:nd)[[ei]]
      dyn.ei = dyn[[paste0("eta_",ei)]]

      res = prod(sapply(1:S.sz,FUN=function(p){
        Yp = ObsModel.transfo$linkS[p]
        sig = Serr[[Yp]]
        tpi = Sobs_i[[Yp]]$time
        Ypi = Sobs_i[[Yp]][,Yp]

        trs = ObsModel.transfo$S[which(ObsModel.transfo$linkS==Yp)]
        Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])

        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
      }))*(1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }else{
    S.margDensity = setNames(sapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
      (1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }

  # s'k(Ri(krij))A(Y|theta0)
  S2.margDensity <- lapply(1:R.sz,FUN=function(k){
    setNames(sapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("eta_",ei)]]

      yGk = ObsModel.transfo$linkR[k]
      sig = Rerr[[yGk]]
      tki = Robs_i[[yGk]]$time
      Zki = Robs_i[[yGk]][,yGk]

      a0 = theta$alpha0[[yGk]]
      a1 = alpha1[[yGk]]
      trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
      Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

      return(1/sig**2*sum((Zki-a0)*Rki)*S.margDensity[[ei]])
    }),paste0("eta_",1:nd))})


  indi.contrib <- sapply(1:R.sz,FUN=function(k){
    sum(mh.parm$Weights*S2.margDensity[[k]])
  })/sum(mh.parm$Weights*S.margDensity)

  return(indi.contrib)
}

# log-lik gradient in alpha1=0 for maximum penalization parameter to reach ----------
lambda.max  <- function(
    dynFUN,
    y,
    mu=NULL,
    Omega=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs=NULL,
    Robs=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    n = NULL,
    prune=NULL,
    parallel = TRUE,
    ncores=NULL,
    onlyLL=FALSE,
    verbose=TRUE){

  if(is.null(data)){
    test <- sapply(c("mu","Omega","theta","alpha1","covariates","ParModel.transfo","ParModel.transfo.inv","Sobs","Robs","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu) || is.null(Omega)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    for(d in 1:length(data)){
      test <- eval(parse(text=(paste0("is.null(",names(data)[d],")"))))
      if(test){
        eval(parse(text=paste0(names(data[d]),"<- data[[d]]")))
      }
    }
  }
  if(is.null(n)){
    n <- floor(100**(1/length(theta$psi_pop)))
  }

  if(parallel){
    if(!is.null(ncores)){
      cluster <- parallel::makeCluster(ncores)
    }else{
      cluster <- parallel::makeCluster(parallel::detectCores())
    }
    doSNOW::registerDoSNOW(cluster)
    }

  N=length(mu)


  i = 1
  if(verbose){
    ntasks <- N
    pb <- utils::txtProgressBar(max = ntasks, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{
    opts <- NULL
  }

  res = foreach::foreach(i = 1:N,.packages = "REMix",.export = c("lambda.max.ind","amgauss.hermite"),.options.snow=opts)%dopar%{
    if(0 %in% diag(Omega[[i]])){
      diag(Omega[[i]])[diag(Omega[[i]])==0] <- 10**(-5)
    }
    lambda.max.ind(mu_i = mu[[i]],
                   Omega_i = Omega[[i]],
                   theta = theta,
                   alpha1 = alpha1,
                   dynFUN = dynFUN,
                   y = y,
                   covariates_i = covariates[i,,drop=F],
                   ParModel.transfo = ParModel.transfo,
                   ParModel.transfo.inv = ParModel.transfo.inv,
                   Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==i,]}),
                   Robs_i = lapply(Robs,FUN=function(R){R[R$id==i,]}),
                   Serr = Serr,
                   Rerr = Rerr,
                   ObsModel.transfo = ObsModel.transfo,
                   n = n,
                   prune = prune)
    }

  if(verbose)
    close(pb)
  if(parallel)
    snow::stopCluster(cluster)
  return(max(Reduce("+",res)))
}

#' Gauss-Hermite approximation of individual log-likelihood derivates
#'
#' @description
#' Compute Gauss-Hermite approximation of individual log-likelihood and its derivates in NLMEM with latent observation process.
#'
#' @details
#' Suppose that we have a differential system of equations containing variables $(S_{p})_{p\leq P}$ and $R$, that depends on some parameters, these dynamics are described by `dynFun`. We write the individual process,  for the individual considered here, over the time, resulting from the differential system for a set of parameters \eqn{\phi_i,(\psi_{li})_{l\leq m}} for the considered individual, as \eqn{S_{p}(\cdot,\phi_i,(\psi_{li})_{l\leq m})=S_{pi}(\cdot)}, \eqn{p\leq P} and \eqn{R(\cdot,\phi_i,(\psi_{li})_{l\leq m})=R_i(\cdot)}. Paremeters are described as \deqn{h_l(\psi_{li}) = h_l(\psi_{lpop})+X_i\beta_l + \eta_{li}} with the covariates of individual \eqn{i}  contains in `covariates_i`, random effects \deqn{\eta_i=(\eta_{li})_{l\leq m}\overset{iid}{\sim}\mathcal N(\mu_i,\Omega_i)} where \eqn{\mu_i} (`mu_i`) is the estimated random effects of individual \eqn{i} and \eqn{\Omega_i} (`Omega_i`) is the diagonal matrix of estimated strandard deviation of random effects of individual \eqn{i}. The population parameters \eqn{\psi_{pop}=(\psi_{lpop})_{l\leq m}} (contain in `theta` through `psi_pop`)  and \eqn{\beta=(\beta_l)_{l\leq m}} (contain in `theta` through `beta`) is the vector of covariates effects on parameters.
#' The rest of the population parameters of the structural model, that hasn't random effetcs, are denoted by \eqn{\phi_i}, and are defined as \eqn{\phi_i=\phi_{pop} + X_i \gamma}, they can depends on covariates effects or be constant over the all population.
#' We assume that individual trajectories \eqn{(S_{pi})_{p\leq P}} are observed through a direct observation model, up to a transformation \eqn{g_p}, \eqn{p\leq P}, at differents times \eqn{(t_{pij})_{p\leq P,j\leq n_{ip}}} : \deqn{Y_{pij}=g_p(S_{pi}(t_{pij}))+\epsilon_{pij}}, contain in `Sobs_i` with error \eqn{\epsilon_p=(\epsilon_{pij})\overset{iid}{\sim}\mathcal N(0,\varsigma_p^2)} for \eqn{p\leq P} (contain in `Serr`).
#' The individual trajectory \eqn{R_{i}} is observed through latent processes, up to a transformation \eqn{s_k}, \eqn{k\leq K}, observed in \eqn{(t_{kij})_{k\leq K,j\leq n_{kij}}} : \deqn{Z_{kij}=\alpha_{k0}+\alpha_{k1} s_k(R_i(t_{kij}))+\varepsilon_{kij}} (contains in arguments `Robs_i`) where \eqn{\varepsilon_k\overset{iid}{\sim} \mathcal N(0,\sigma_k^2)} (contains in `Rerr`).
#'
#' @param mu_i Named vector of individual random effetcs estimation (size m);
#' @param Omega_i Named diagonal matrix of individual standard deviation estimation (size m x m) ;
#' @param theta model parameter ; contain phi_pop (size L), psi_pop (size m), gamma (list of size L, containing vector of size m), beta (list of size m, containing vector of size m), alpha0 (size K), and omega (size m) ;
#' @param alpha1 Regularization parameter (size K ) ; must be in the same order as in Robs_i.
#' @param dynFun Dynamic function ;
#' @param y Initial condition of the model, conform to what is asked in dynFun ;
#' @param covariates_i Covariates line matrix of the individual ;
#' @param ParModel.transfo named list of transformation function for individual parameter model, (the name must be consistent with phi_pop) (size <=L, missing is set to identity) ;
#' @param ParModel.transfo.inv named list of inverse transformation function for individual parameter model, (the name must be consistent with phi_pop) (size <=L, missing is set to identity) ;
#' @param Sobs_i list of direct observation (size P), each element contain time and observation (in column 3) ;
#' @param Robs_i list of latent observation (size K), each element contain time and observation (in column 3) ;
#' @param Serr vector of size P containing estimated error model constant (must be in the same order os in Sobs_i)
#' @param Rerr vector of size K containing estimated error model constant (must be in the same order os in Robs_i)
#' @param ObsModel.transfo
#' list of 2 list of P,K transformation (need to include identity transformation), named with `S` and `R` :
#'
#'   - ObsModel.transfo$S correspond to the transformation used for direct observation model. For each \eqn{Y_p=h_p(S_p)} the order (as in `Sobs_i`) must be respected and the name indicated which dynamic from dynFun is observed through this variables \eqn{Y_p};
#'
#'   - ObsModel.transfo$R correspond to the transformation used for the latent process, as it is now, we only have one latent dynamic so necessarily \eqn{s_k} is applied to `R` but for each \eqn{Y_k} observed, transformation could be different so need to precise as many as in `Robs_i` ; the name need be set to precise the dynamic from dynFun to identify the output.
#' @param n (default floor(100**(1/length(theta$psi_pop))) number of points for gaussian quadrature ;
#' @param prune (default NULL) percentage in [0;1] for prunning.
#' @param onlyLL (default FALSE) if only the log-likelihood need to be return.
#'
#' @return a list with the approximation by Gauss-Hermite quadrature of the likelihood `Li`, the log-likelihoo `LLi`, the gradient of the log-likelihood `dLLi` and the Hessien of the log-likelihood `ddLLi` in point \eqn{\theta,\alpha}.
#' @export
#' @seealso [gh.int()]
#' @import stats
#'
#' @examples
#' #TODO

gh.LL.ind <- function(
    dynFun,
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
  # check to do

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
  # N = nrow(covariates_i)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

  # Need to compute, for each eta_i x individual i, the dynamics of the model
  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
    dyn_eta_i <- dynFun(all.tobs,y,unlist(unname(PSI_i)))

    return(dyn_eta_i)
  }),paste0("eta_",1:nd))

  # Compute marginal latent density for each individual and value of RE
  R.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("eta_",ei)]]


      res = prod(sapply(1:R.sz,FUN=function(k){
        sig = Rerr[[k]]
        tki = Robs_i[[k]]$time
        Zki = Robs_i[[k]][,3]

        a0 = theta$alpha0[[k]]
        a1 = alpha1[[k]]
        trs = ObsModel.transfo[[2]][k]
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
        sig = Serr[[p]]
        tpi = Sobs_i[[p]]$time
        Ypi = Sobs_i[[p]][,3]

        trs = ObsModel.transfo[[1]][p]
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
        sig = Rerr[[k]]
        tki = Robs_i[[k]]$time
        Zki = Robs_i[[k]][,3]

        a0 = theta$alpha0[[k]]
        a1 = alpha1[[k]]
        trs = ObsModel.transfo[[2]][k]
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
        sig = Rerr[[k]]
        tki = Robs_i[[k]]$time
        Zki = Robs_i[[k]][,3]

        a0 = theta$alpha0[[k]]
        a1 = alpha1[[k]]
        trs = ObsModel.transfo[[2]][k]
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

              sig = Rerr[[k]]
              tki = Robs_i[[k]]$time
              Zki = Robs_i[[k]][,3]

              a0 = theta$alpha0[[k]]
              a1 = alpha1[[k]]
              trs = ObsModel.transfo[[2]][k]
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

  ## Transformation pour N(mu,sigma)  (sigma=sd)
  ww <- gh$Weights * sqrt(2) * sd * exp(gh$Points**2)
  xx <- mu + sd * sqrt(2) * gh$Points

  return(data.frame(Points=xx,Weights=ww))
}


# Multi-dimensional Adaptative Hermite Gauss -----------------------------------------
amgauss.hermite <- function(n,mu=rep(0,ncol(Omega)),Omega=diag(rep(1,length(mu))),prune=NULL){
  # L'idée est on prends les points en dim 1 et on les rotatent/translatent
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


################# lambda_max ##################################
lambda.max.ind <- function(
    dynFun,
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
  # N = nrow(covariates_i)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

  # Need to compute, for each eta_i x individual i, the dynamics of the model
  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
    dyn_eta_i <- dynFun(all.tobs,y,unlist(unname(PSI_i)))

    return(dyn_eta_i)
  }),paste0("eta_",1:nd))

  # Compute marginal latent density for each individual and value of RE
  R.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
    dyn.ei = dyn[[paste0("eta_",ei)]]


    res = prod(sapply(1:R.sz,FUN=function(k){
      sig = Rerr[[k]]
      tki = Robs_i[[k]]$time
      Zki = Robs_i[[k]][,3]

      a0 = theta$alpha0[[k]]
      a1 = alpha1[[k]]
      trs = ObsModel.transfo[[2]][k]
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
        sig = Serr[[p]]
        tpi = Sobs_i[[p]]$time
        Ypi = Sobs_i[[p]][,3]

        trs = ObsModel.transfo[[1]][p]
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

      sig = Rerr[[k]]
      tki = Robs_i[[k]]$time
      Zki = Robs_i[[k]][,3]

      a0 = theta$alpha0[[k]]
      a1 = alpha1[[k]]
      trs = ObsModel.transfo[[2]][k]
      Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

      return(1/sig**2*sum((Zki-a0)*Rki)*S.margDensity[[ei]])
  }),paste0("eta_",1:nd))})


  indi.contrib <- sapply(1:R.sz,FUN=function(k){
    sum(mh.parm$Weights*S2.margDensity[[k]])
  })/sum(mh.parm$Weights*S.margDensity)

return(indi.contrib)
}

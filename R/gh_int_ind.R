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
#' @param theta model parameter ; contain phi_pop (size L), psi_pop (size m), gamma (list of size L, containing vector of size m), beta (list of size m, containing vector of size m), alpha0 (size K);
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
#' @param ObsModel.transfo list of 2 list of P,K transformation (need to include identity transformation), named with `S` and `R` :
#'
#'   - ObsModel.transfo$S correspond to the transformation used for direct observation model. For each \eqn{Y_p=h_p(S_p)} the order (as in `Sobs_i`) must be respected and the name indicated which dynamic from dynFun is observed through this variables \eqn{Y_p};
#'
#'   - ObsModel.transfo$R correspond to the transformation used for the latent process, as it is now, we only have one latent dynamic so necessarily \eqn{s_k} is applied to `R` but for each \eqn{Y_k} observed, transformation could be different so need to precise as many as in `Robs_i` ; the name need be set to precise the dynamic from dynFun to identify the output.
#' @param n (default floor(100**(1/length(theta$psi_pop))) number of points for gaussian quadrature ;
#' @param prune (default NULL) percentage in [0;1] for prunning.
#'
#' @return a list with the approximation by Gauss-Hermite quadrature of the likelihood `Li`, the log-likelihoo `LLi`, the gradient of the log-likelihood `dLLi` and the Hessien of the log-likelihood `ddLLi` in point \eqn{\theta,\alpha}.
#' @export
#' @seealso [gh.int()]
#' @import stats
#'
#' @examples
#' #TODO

gh.int.ind <- function(mu_i,Omega_i,theta,alpha1,dynFun,y,covariates_i,ParModel.transfo,ParModel.transfo.inv,Sobs_i,Robs_i,Serr,Rerr,ObsModel.transfo,n=floor(100**(1/length(theta$psi_pop))),prune=NULL){
  # check to do

  if(ncol(Omega_i)!=1){
    mh.parm <- mgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
  }else{
    mh.parm <- gauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
  }
  nd = length(mh.parm$Weights)
  N = nrow(covariates_i)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  root.Omega = solve(Omega_i)
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
      dyn.ei = dyn[[paste0("eta_",ei)]]

      res = prod(sapply(1:S.sz,FUN=function(p){
        sig = Serr[[p]]
        tpi = Sobs_i[[p]]$time
        Ypi = Sobs_i[[p]][,3]

        trs = ObsModel.transfo[[1]][p]
        Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])

        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
      }))
    }),paste0("eta_",1:nd))
  }else{
    S.margDensity = setNames(rep(1,nd),paste0("eta_",1:nd))
  }


  # Compute individual log-Likelihood
  Li = sum(mh.parm$Weights*R.margDensity*S.margDensity)
  LLi = log(Li)

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

      return(1/(sig**2)*(sum(alpha1[[k]]*Rki*(Zki-theta$alpha0[[k]]-alpha1[[k]]*Rki))))
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
        (1/sig**4*(sum(alpha1[[k]]*Rki*(Zki-theta$alpha0[[k]]-alpha1[[k]]*Rki)))**2 +
          1/sig**2*(sum(Rki*(Zki-theta$alpha0[[k]]-2*alpha1[[k]]*Rki))))
        )
    }),paste0("eta_",1:nd))
  })

  diag(ddLLi) = sapply(ddfact,FUN=function(f){
    return(1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
  })-dLLi**2

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
              1/sig**2*(sum(alpha1[[k]]*Rki*(Zki-theta$alpha0[[k]]-alpha1[[k]]*Rki)))**2
            )
        }))
        }),paste0("eta_",1:nd))

        ddLLi[k1,k2] <- ddLLi[k2,k1] <- 1/Li*sum(mh.parm$Weights*ddfact*R.margDensity*S.margDensity)-prod(dLLi[c(k1,k2)])
      }
    }
  }

  return(list(Li=Li,LLi=LLi,dLLi=dLLi,ddLLi=ddLLi))
}


# Hermite polynomial ------------------------------------------------------
# \int f(x)exp(-x²)dx=sum wi * f(xi)
hermite <- function(n){
  # Hn+1 = XHn - nHn-1
  # degré max n
  # n>=1
  H0 <- c(1,rep(0,n))
  H1 <- c(0,2,rep(0,n-1))
  if(n>1){
    for(k in 1:(n-1)){
      H0.prev = H0
      H0 <- H1
      H1 <- 2*(data.table::shift(H1,fill = 0) - k*H0.prev)
      # print(H1)
    }
  }
  # H0 c'est notre Hn-1 et H1 notre H1

  ## Polynôme de Gauss
  x <- sort(Re(polyroot(H1))) # calcul des xi comme racine du polynôme d'hermite Hn(x)
  H0x <- sapply(x,FUN=function(xi){
    res <- 0
    for(k in 0:(n-1)){
      res <- res + H0[k+1]*xi**k
    }
    return(res)
  }) # calcul du polynôme Hn-1 en chaque points (xi)
  w <- 2**(n + 1) * factorial(n) * sqrt(pi)/(2*n*H0x)^2 # calcul des poids wi= 2**(n+1)*factorial(n)*sqrt(pi) / H'n(xi)**2 et H'n(x)=2nHn-1(x)

  return(data.frame(Points=x,Weights=w))
}
# tester pour n=10 ; f=id ; f=x ; f=x^2 et à l'air good


# One dimensional hermite gauss  ------------------------------------------
# changement de variable pour int f(eta) deta selon mesure normale i.e int f(x) 1/(sd**2 sqrt(2*pi))*exp(-1/2((x-mu)/sd)**2)dx = sum(wi f(xi))
gauss.hermite <- function(n,mu=0,sd=1){
  gh <- hermite(n)

  ## Transformation pour N(mu,sigma)  (sigma=sd)
  ww <- gh$Weights/sqrt(pi)
  xx <- mu + sd * sqrt(2) * gh$Points

  return(data.frame(Points=xx,Weights=ww))
}
# tester pour n=10 ; f=id ; f=x ; f=x^2 et à l'air good


# Multi-dimensional Hermite Gauss -----------------------------------------
# ref http://www.jaeckel.org/ANoteOnMultivariateGaussHermiteQuadrature.pdf
# The notion of ‘pruning’ can be used to eliminate those points with very small weights.
mgauss.hermite <- function(n,mu=rep(0,ncol(Omega)),Omega=diag(rep(1,length(mu))),prune=NULL){
  # L'idée est on prends les points en dim 1 et on les rotatent/translatent
  dm = length(mu)

  gh  <- gauss.hermite(n)
  # plot(idx)  dm  <- length(mu)
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
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  # colnames(pts) <-c("X","Y")
  # ggplot(mapping=aes(x=X,y=Y))+geom_density_2d(data=data,aes(color = after_stat(level)))+scale_color_viridis_c()+geom_point(data=as.data.frame(pts,wts=wts),aes(size=-5/log(wts)))
  return(list(Points=pts, Weights=wts))
}
# j'ai bien la densité d'une loi normale vaut 1
# j'arrive à retrouver P(X<0,Y<0) pour loi normale centrée réduite
# idem en plus grande dimension n=6 au lieu de 2
# > Omega=diag(rep(1,6))
# > n
# [1] 2
# > mu=rep(0,6)
# > t <- mgauss.hermite(n,mu,Omega)
# > sum(t$Weights*apply(t$points,1,FUN=gfun))
# [1] 0.015625
# > 1/(2**6)
# [1] 0.015625

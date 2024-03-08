#' Title
#'
#' @param Omega a
#' @param theta a
#' @param alpha1 a
#' @param dynFun a
#' @param y a
#' @param covariates a
#' @param ParModel.transfo a
#' @param ParModel.transfo.inv aa
#' @param Sobs a
#' @param Robs a
#' @param Serr a
#' @param Rerr a
#' @param ObsModel.transfo a
#' @param ObsModel.transfo.inv a
#' @param n a
#' @param prune a
#'
#' @return a
#' @import foreach
#' @export
#'
#' @examples a
gh.int.ind <- function( Omega, # NAMED matrix of covariance, size m
                        theta, # model parameter contain Psi_pop (size L), psi_pop (size m), Beta (size L), beta size(m), alpha0 (size K)
                        alpha1, # regularization parameter (size K )
                        dynFun, # Function dynamic
                        y,      # initial condition, confom to what is asked in dynFun
                        covariates, # covariance matrix, individual in line (N x n )
                        ParModel.transfo, # list of transfo function, named list (consistent with Psi_pop) size <= L
                        ParModel.transfo.inv, # list of transfo function, named list (consistent with psi_pop), size <= K
                        Sobs, # list (size P) with matrix of time and all observation in column (N+1 column)
                        Robs, # list (size K) with matrix of time and all observation in column (N+1 column) ; with NA when individual not observed (?)
                        Serr, # vector (size P) of sigma
                        Rerr, # vector (size K) of sigma
                        ObsModel.transfo, # list of 2 list of P,K transfo (need to include id) NAMED
                        # ObsModel.transfo$S correspond to the transfo used for direct observation model. FOr each Yp=hp(Sp) the order (as in Sobs) must be respected and the name indicated which dynamic for dynFun is observed through this variables Yp
                        # ObsModel.transfo$R correspond to the transfo used for the latent process, as it is now, we only have on latent dynamic so necessarily sk is applied to R but for each Yk observed, trans could be different so need to precise as many as in Robs ; need to precise R dyn name to identify it in the output of dynFun
                        # DO NOT FORGET TO INCLUDE THE IDENTITY TRANSFORMATION
                        n=3,
                        prune=NULL){
  # doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()-5))
  # check to do

  if(ncol(Omega)!=1){
    mh.parm <- mgauss.hermite(n,Omega=Omega,prune=prune)
  }else{
    mh.parm <- gauss.hermite(n)
  }
  nd = nrow(mh.parm$points)
  N = nrow(covariates)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  root.Omega = solve(Omega)
  all.tobs = sort(union(unlist(lapply(Robs,FUN=function(x){x$time})),unlist(lapply(Sobs,FUN=function(x){x$time}))))

  # Need to compute, for each eta_i x individual i, the dynamics of the model
  dyn <- setNames(foreach(i = 1:N,.combine=append,.packages="LassoLLPen") %dopar% {
    list(setNames(lapply(split(mh.parm$points,1:nd),FUN=function(eta_i){
      PSI_i  = indParm(theta[c("Psi_pop","psi_pop","Beta","beta")],covariates[i,,drop=F],setNames(eta_i,colnames(Omega)),ParModel.transfo,ParModel.transfo.inv)
      dyn_eta_i <- dynFun(all.tobs,y,unlist(unname(PSI_i)))

      return(dyn_eta_i)
    }),paste0("eta_",1:nd)))
  }, paste0("ind_",1:N))

  # Compute marginal latent density for each individual and value of RE
  R.margDensity <- foreach(i = 1:N,.combine = rbind,.packages = "LassoLLPen") %dopar%{
    # For each individual and eta_i fixed, we want to compute the marginal density of (R_k)_{k\leq K}
    eval_eta_i = setNames(sapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("ind_",i)]][[paste0("eta_",ei)]]

      res = prod(sapply(1:R.sz,FUN=function(k){
        sig = Rerr[[k]]
        Zki = Robs[[k]][,c("time",paste0("ind_",i))]
        tki = Zki[!is.na(Zki[,2]),1]
        Zki = Zki[Zki$time %in% tki,2]

        a0 = theta$alpha0[[k]]
        a1 = alpha1[[k]]
        trs = ObsModel.transfo[[2]][k]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])
        plot(Rki)
        plot(Zki)

        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2))
      }))
    }),paste0("eta_",1:nd))
  }

  # Compute marginal density for each individual and value of RE
  S.margDensity <- foreach(i = 1:N,.combine = rbind,.packages = "LassoLLPen") %dopar%{
    # For each individual and eta_i fixed, we want to compute the marginal density of (S_p)_{p\leq P}
    eval_eta_i = setNames(sapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("ind_",i)]][[paste0("eta_",ei)]]

      res = prod(sapply(1:S.sz,FUN=function(p){
        sig = Serr[[p]]
        Ypi = Sobs[[p]][,c("time",paste0("ind_",i))]
        tpi = Ypi[!is.na(Ypi[,2]),1]
        Ypi = Ypi[Ypi$time %in% tpi,2]

        trs = ObsModel.transfo[[1]][p]
        Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])
        plot(Ypi)
        plot(Spi)

        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
      }))
    }),paste0("eta_",1:nd))
  }

  # Compute individual Likelihood
  LLi = sapply(1:N,FUN=function(i){
    sum(mh.parm$weights*R.margDensity[i,]*S.margDensity[i,])
  })

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
# > sum(t$weights*apply(t$points,1,FUN=gfun))
# [1] 0.015625
# > 1/(2**6)
# [1] 0.015625

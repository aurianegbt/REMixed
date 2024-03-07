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
#' @export
#'
#' @examples a
gh.int.ind <- function( Omega, # NAMED matrix of covariance, size m
                        theta, # model parameter contain Psi_pop (size L), psi_pop (size m), Beta (size L), beta size(m), alpha0 (size K)
                        alpha1, # regularization parameter (size K )
                        dynFun, # Function dynamic
                        identifier.SR, # SAME ORDER AS IN OBSMODEL[[1]] and Sobs and Serr
                        y,      # initial condition, confom to what is asked in dynFun
                        covariates, # covariance matrix, individual in line (N x n )
                        ParModel.transfo, # list of transfo function, named list (consistent with Psi_pop) size <= L
                        ParModel.transfo.inv, # list of transfo function, named list (consistent with psi_pop), size <= K
                        Sobs, # list (size P) with matrix of time and all observation in column (N+1 column)
                        Robs, # list (size K) with matrix of time and all observation in column (N+1 column) ; with NA when individual not observed (?)
                        Serr, # vector (size K) of sigma
                        Rerr, # vector (size P) of sigma
                        ObsModel.transfo, # list of 2 list of P,K transfo (need to include id)
                        n=3,
                        prune=NULL){
  # doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()-5))
  # check to do

  if(ncol(Omega)!=1){
    mh.parm <- mgauss.hermite(n,Omega,prune=prune)
  }else{
    mh.parm <- gauss.hermite(n)
  }
  nd = nrow(mh.parm$points)
  N = nrow(covariates)
  root.Omega = solve(Omega)
  all.tobs = sort(union(unlist(lapply(Robs,FUN=function(x){x$time})),unlist(lapply(Sobs,FUN=function(x){x$time}))))

  # Need to compute, for each eta_i x individual i, the dynamics of the model
  dyn <- foreach(i = 1:N,.combine=append,.packages="LassoLLPen") %dopar% {
    list(lapply(split(mh.parm$points,1:nd),FUN=function(eta_i){
      PSI_i  = indParm(theta[c("Psi_pop","psi_pop","Beta","beta")],covariates[i,,drop=F],setNames(eta_i,colnames(Omega)),ParModel.transfo,ParModel.transfo.inv)
      dyn_eta_i <- dynFun(all.tobs,y,unlist(unname(PSI_i)))

      return(dyn_eta_i)
    }))
  }
  names(dyn) <- paste0("ind_",1:N)

  ## Calcul de f(Zi;theta_l,alpha_l,eta_i) et A(Yi;theta_l,eta_i)
  # eta_i en colonne, dans l'order, ind en ligne
  f.compute_i <- foreach(i = 1:N,.combine=rbind,.packages="LassoLLPen") %dopar% {
    sapply(1:nd,FUN=function(ei){
      Zi = lapply(Robs,FUN=function(x){x[,c(1,i+1)]})
      dyn_ei = dyn[[i]][[ei]][,c("time",identifier.SR$R)]

      aux = prod(sapply(1:length(Zi),FUN=function(k){
        tkij = Zi[[k]][!is.na(Zi[[k]][[2]]),"time"]
        Zki = Zi[[k]][!is.na(Zi[[k]][[2]]),2]
        R.tkij = dyn_ei[dyn_ei[,1] %in% tkij,2]

        return(prod(
          1/(Rerr[[k]]*sqrt(2*pi))*
            exp(-1/2*((Zki-theta$alpha0[k]-alpha1[k]*ObsModel.transfo[[2]][[k]](R.tkij))/Rerr[[k]])**2)))
      }))
      return(as.numeric(aux)) # aux is f(Zi;theta_l,alpha_l,eta_i) for a i and a eta_i
    })
  }

  A.compute_i <- foreach(i = 1:N,.combine=rbind,.packages="LassoLLPen") %dopar% {
    sapply(1:nd,FUN=function(ei){
     Yi = lapply(Sobs,FUN=function(x){x[,c(1,i+1)]})
     dyn_ei = dyn[[i]][[ei]][,c("time",identifier.SR$S)]

      aux = prod(sapply(1:length(Yi),FUN=function(p){
        tpij = Yi[[p]][!is.na(Yi[[p]][[2]]),"time"]
        Ypi = Yi[[p]][!is.na(Yi[[p]][[2]]),2]
        S.tpij = dyn_ei[dyn_ei[,1] %in% tpij,p+1] # identifier.SR same order as in Yp

        return(prod(
          1/(Serr[[p]]*sqrt(2*pi))*exp(-1/2*((Ypi-ObsModel.transfo[[1]][[p]](S.tpij))/Serr[[p]])**2)))})) * 1/((2*pi)**(ncol(Omega)/2)*det(Omega)**(1/2))*exp(-1/2*matrix(mh.parm$points[ei,],nrow=1)%*%root.Omega%*%matrix(mh.parm$points[ei,],ncol=1))
      return(as.numeric(aux)) # aux is A(Yi;theta_l,eta_i) for a i and a eta_i
    })
  }
  ## Calcul de LL(theta_l,alpha_l)
  LLi = sapply(1:N,FUN=function(i){
    int_i <- sum(mh.parm$weights * f.compute_i[i,]*A.compute_i[i,])
    # int_i <- sum(mh.parm$weights * mpfr(f.compute_i[i,],precBits = 100)*mpfr(A.compute_i[i,],precBits = 100))
  })
  LL = log(prod(LLi)) # Il y a des -inf... comment les gérer ?

  ## Calcul du gradient  dLL(theta_l,alpha_l)
  dLL = sapply(1:length(alpha1),FUN=function(k){





    return(dLL_k)
  })

  ## Calcul de d²LL(theta_l,alpha_l)

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

set.seed(1710)
N=100

# R.E. argument -----------------------------------------------------------
omega_fM2 = 0.8
omega_theta = 0.5
Omega = diag(c(omega_fM2**2,omega_theta**2)) # m=2
rownames(Omega) <- colnames(Omega) <- c("fM2","theta")
# Parameter value  --------------------------------------------------------
Psi_pop = c(delta_V=2.7,delta_S=0.01,delta_Ab=0.03,t_0=0,t_inj=21)
psi_pop = c(fM2 = 4.5,theta=18.7) # m=2
Beta = NULL
beta = NULL
alpha0 = rep(0,2) # K=2

theta=setNames(list(Psi_pop,psi_pop,Beta,beta,alpha0),c("Psi_pop","psi_pop","Beta","beta","alpha0"))

alpha1 = c(1,1) # K=2
# Dynamic function  -------------------------------------------------------
dynFun = Clairon
identifier.SR = list(S=c("A"),R=c("S")) # S for state, R for i don't know the f**k its name should be ; order must be the same as in Sobs !!!!
y = c(S=0,A=0.5)
covariates =mvtnorm::rmvnorm(N,mean = rep(0,10)) # N=100
ParModel.transfo = list(fM2=log,theta=log) # m=2
ParModel.transfo.inv = list(fM2 = exp,theta=exp) # m=2
# Observation -------------------------------------------------------------
Sobs = data.frame(time = c(0,2,7,14,21,28,35)) # P = 1
Robs = lapply(1:2,FUN=function(x){data.frame(time = c(1:14,21:35))}) # K=2
random.effect <- mvtnorm::rmvnorm(N,sigma=Omega)


for(i in 1:N){
  par.i <- indParm(theta,covariates,setNames(random.effect[i,],c("fM2","theta")),ParModel.transfo,ParModel.transfo.inv)

  dyn <- dynFun(t=sort(union(c(0,2,7,14,21,28,35),c(1:14,21:35))),y,unlist(unname(par.i)))

  Sobs[,paste0("ind_",i)] <- log10(dyn[dyn[,1] %in% Sobs$time,"A"])+rnorm(length(Sobs$time),sd=0.23)
  for(k in 1:2){
    Robs[[k]][,paste0("ind_",i)] <- theta$alpha0[k] + alpha1[k]*dyn[dyn[,1] %in% Robs[[k]]$time,"S"] + rnorm(length(Robs[[k]]$time),sd=0.1)
  }
}
Sobs <- list(Sobs)

Serr = 0.23 # P=1 sd here !!!
Rerr = c(0.1,0.1) #K=2 sd here !!!
ObsModel.transfo = list(S=list(A=log10),R=list(S=function(x){x},S=function(x){x})) # P+K=3

# Gaussian quadrature parameter -------------------------------------------
n = 5
prune=NULL


# Add non observed --------------------------------------------------------

pb_missing = 0.01
# Robs missing :
for(k in 1:length(Robs)){
  miss = lapply(1:N,FUN=function(x){rbinom(length(Robs[[k]]$time),size = 1,prob=pb_missing)})
  for( i in 1:N){
    Robs[[k]][,i+1][as.logical(miss[[i]])] <- NA
  }
}

# Sobs missing :
for(p in 1:length(Sobs)){
  miss = lapply(1:N,FUN=function(x){rbinom(length(Sobs[[p]]$time),size = 1,prob=pb_missing)})
  for( i in 1:N){
    Sobs[[p]][,i+1][as.logical(miss[[i]])]  <- NA
  }
}

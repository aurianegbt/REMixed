set.seed(1710)
N=100

# R.E. argument -----------------------------------------------------------
omega_ka = 0.8
Omega = matrix(c(omega_ka**2)) # m=2
rownames(Omega) <- colnames(Omega) <- c("ka")
# Parameter value  --------------------------------------------------------
Psi_pop = NULL
psi_pop = c(ka=0.7) # m=2
Beta = NULL
beta = NULL
alpha0 = rep(0,2) # K=2

theta=setNames(list(Psi_pop,psi_pop,Beta,beta,alpha0),c("Psi_pop","psi_pop","Beta","beta","alpha0"))

alpha1 = c(1,2) # K=2
# Dynamic function  -------------------------------------------------------
dynFun = PK
identifier.SR = list(S=c(),R=c("C")) # S for state, R for i don't know the f**k its name should be ; order must be the same as in Sobs !!!!
y = c(C=100)
covariates =mvtnorm::rmvnorm(N,mean = rep(0,10)) # N=100
ParModel.transfo = list(ka=log) # m=2
ParModel.transfo.inv = list(ka=exp) # m=2
# Observation -------------------------------------------------------------
Sobs = data.frame(time = c()) # P = 0
Robs = lapply(1:2,FUN=function(x){data.frame(time = c(0,0.5,1,2,5,10,20,50,100,150))}) # K=2
random.effect <- mvtnorm::rmvnorm(N,sigma=Omega)

for(i in 1:N){
  par.i <- indParm(theta,covariates,setNames(random.effect[i,],c("ka")),ParModel.transfo,ParModel.transfo.inv)

  dyn <- dynFun(t=c(0,0.5,1,2,5,10,20,50,100,150),y,unlist(unname(par.i)))

  # Sobs[,paste0("ind_",i)] <- c()
  for(k in 1:2){
    Robs[[k]][,paste0("ind_",i)] <- theta$alpha0[k] + alpha1[k]*dyn[dyn[,1] %in% Robs[[k]]$time,"C"] + rnorm(length(Robs[[k]]$time),sd=0.1)
  }
}
Sobs <- list(Sobs)

Serr = NULL # P=1 sd here !!!
Rerr = c(0.1,0.1) #K=2 sd here !!!
ObsModel.transfo = list(list(),list(function(x){x},function(x){x})) # P+K=3
ObsModel.transfo.inv = list(list(),list(function(x){x},function(x){x})) # P+K=3

# eta_i =
# [1]  0.054996509 -0.034226751 -0.477703301 -0.451741725 -0.043533260 -0.226711587
# [7] -0.398359707  2.877128236 -0.258358112 -0.254430137  0.878445907 -0.555149047
# [13] -0.342787977  2.213159979 -0.464823397 -0.524048880  1.689045684 -0.558138743
# [19] -0.144577618  0.180446484 -0.443648334  0.925238797  2.835551599  0.184879168
# [25] -0.228151684  0.104802398 -0.444859995  0.063846276 -0.307536977  0.544182144
# [31]  0.244669675  0.426791513  0.003417717 -0.443268425 -0.417495368 -0.392370759
# [37]  0.077695584  0.590813250  0.770689435  1.088752679  0.886811716  0.113178865
# [43] -0.210966791  0.100587052  0.378789267  1.365335927 -0.160557684 -0.083353777
# [49] -0.252353390 -0.351591141 -0.081257700  0.010349836 -0.144239431 -0.512953110
# [55] -0.082299604 -0.358291066 -0.418587540 -0.455704131  0.869646011 -0.642447449
# [61] -0.358571412  0.890211015 -0.417399946  0.261452277  0.537897901 -0.558512177
# [67]  1.874420484 -0.249705543 -0.136791178 -0.161171567 -0.547270240 -0.547678554
# [73] -0.516684867  0.188947835 -0.603483278  0.294252130  0.477078892  0.860666361
# [79]  1.911094889  1.548293816  3.169860075 -0.331887277 -0.040131103  0.065730838
# [85] -0.451290933 -0.491747435 -0.415401466  0.225068725 -0.605559423  0.054701343
# [91] -0.234096968  0.881758584  0.502275489 -0.610941449  1.257930535  0.340829044
# [97]  0.435917588  1.369534672  0.168607886 -0.227825424


# Gaussian quadrature parameter -------------------------------------------
n = 100
prune=NULL

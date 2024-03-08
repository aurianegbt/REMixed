Psi_pop = c(delta_S = 0.231, delta_L = 0.000316)
psi_pop = c(delta_Ab = 0.025,phi_S = 3057, phi_L = 16.6)
Beta = NULL

covariates = data.frame(cAGE = runif(1,-15,15), G1 = rnorm(1), G2 = rnorm(1))

beta = list(delta_Ab=c(0,1.2,0),phi_S = c(0.93,0,0),phi_L=c(0,0,0.8))

theta=setNames(list(Psi_pop,psi_pop,Beta,beta),c("Psi_pop","psi_pop","Beta","beta"))

eta_i = c(delta_Ab = rnorm(1,0,0.3),phi_S=rnorm(1,0,0.92),phi_L=rnorm(1,0,0.85))

transfo = list(delta_Ab=log,phi_S=log,phi_L=log)
transfo.inv = list(delta_Ab = exp,phi_S=exp,phi_L=exp)

save(theta,covariates,eta_i,transfo,transfo.inv,file="Exemple/indParmPasin.RData")


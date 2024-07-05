#' Generate trajectory of the Humoral Immune Response to a Prime-Boost Ebola Vaccine
#'
#'The model correspond to the dynamics of the humoral response, from 7 days after the boost immunization. Model is defined as
#' \deqn{ \left\{\begin{matrix}\frac{d}{dt} Ab(t) &=& \theta_S S(t) + \theta_L L(t) - \delta_{Ab} Ab(t)  \\ \frac{d}{dt}S(t) &=& -\delta_S S(t) \\ \frac{d}{dt} L(t) &=& -\delta_L L(t)\end{matrix}\right. }
#'
#' @param t vector of time ;
#' @param y initial condition, named vector of form c(A=A0,S=S0,L=L0) ;
#' @param parms named vector of model parameter ; should contain ""theta_S","theta_L","delta_Ab","delta_S","delta_L".
#'
#' @return Matrix of time and observation of antibody titer Ab, and ASCs S and L.
#' @export
#' @seealso [indParm()]
#'
#' @examples
# y = c(A=0,S=5,L=5)
#
# parms = c(theta_S = 611,
#           theta_L = 3.5,
#           delta_Ab = 0.025,
#           delta_S = 0.231,
#           delta_L = 0.000152)
#
# t = seq(0,100,5)
#
# res <- Pasin(t,y,parms)
#
# plot(res)
#'
#' @references Pasin C, Balelli I, Van Effelterre T, Bockstal V, Solforosi L, Prague M, Douoguih M, Thiébaut R, for the EBOVAC1 Consortium. 2019. Dynamics of the humoral immune response to a prime-boost Ebola vaccine: quantification and sources of variation. J Virol 93:e00579-19. https://doi.org/10.1128/JVI.00579-19.

Pasin <- function(t,y,parms){
  #phi_S,phi_L,delta_Ab,delta_S=0.23,delta_L=0.000316
  if(!setequal(names(y),c("A","S","L"))){
    bonus = setdiff(names(y),c("A","L","S"))
    malus = setdiff(c("A","S","L"),names(y))
    if(length(bonus)!=0 & length(malus) !=0){
    stop(paste0("Missing initial condition for ",malus," and ",bonus," isn't in the model."))
    }else if(length(bonus)!=0){
      stop(paste0(bonus," isn't dynamic of the model."))
    }else if(length(minus) !=0){
      stop(paste0("Missing initial condition for ",malus))
    }
  }
  if(!setequal(names(parms),c("theta_S","theta_L","delta_Ab","delta_S","delta_L"))){
    stop(paste0("Missing parmeters ",setdiff(c("theta_S","theta_L","delta_Ab","delta_S","delta_L"),names(parms))," and ",setdiff(names(parms),c("theta_S","theta_L","delta_Ab","delta_S","delta_L"))," isn't in the model."))
  }
  out <- deSolve::ode(y,t,modelPasin,parms)
  return(out)
}


# Model -------------------------------------------------------------------

modelPasin <- function(t,y,parms){
  with(as.list(c(y, parms)), {
    dA = theta_S*S + theta_L-L-delta_Ab*A
    dS = -delta_S *S
    dL = -delta_L * L

    list(c(dA,dS,dL))
  })
}

#' Generate trajectory of the Humoral Immune Response to a Prime-Boost Ebola Vaccine
#'
#'The model correspond to the dynamics of the humoral response, from 7 days after the boost immunization. Model is defined as
#' \deqn{ \frac{d}{dt} Ab(t) = \varphi_S e^{-\delta_S t} + \varphi_L e^{-\delta_L t} - \delta_{Ab} Ab(t) }
#'
#' @param t vector of time ;
#' @param y initial condition, named vector of form c(A=A0) ;
#' @param parms named vector of model parameter ; should contain ""phi_S","phi_L","delta_Ab","delta_S","delta_L".
#'
#' @return Matrix of time and observation of antibody titer Ab.
#' @export
#'
#' @examples
#' y = c(A=0)
#'
#' parms = c(phi_S = 3057,
#'           phi_L = 16.6,
#'           delta_Ab = 0.025,
#'           delta_S = 0.231,
#'           delta_L = 0.000316)
#'
#' t = seq(0,100,5)
#'
#' res <- Pasin(t,y,parms)
#'
#' plot(res)
#'
#' @references Pasin C, Balelli I, Van Effelterre T, Bockstal V, Solforosi L, Prague M, Douoguih M, Thiébaut R, for the EBOVAC1 Consortium. 2019. Dynamics of the humoral immune response to a prime-boost Ebola vaccine: quantification and sources of variation. J Virol 93:e00579-19. https://doi.org/10.1128/JVI.00579-19.

Pasin <- function(t,y,parms){
  #phi_S,phi_L,delta_Ab,delta_S=0.23,delta_L=0.000316
  if(!identical(names(y),"A")){
    stop(paste0("Missing initial condition for ",setdiff(c("A"),names(y))," and ",setdiff(names(y),c("A"))," isn't in the model."))
  }
  if(!setequal(names(parms),c("phi_S","phi_L","delta_Ab","delta_S","delta_L"))){
    stop(paste0("Missing parmeters ",setdiff(c("phi_S","phi_L","delta_Ab","delta_S","delta_L"),names(parms))," and ",setdiff(names(parms),c("phi_S","phi_L","delta_Ab","delta_S","delta_L"))," isn't in the model."))
  }
  out <- deSolve::ode(y,t,modelPasin,parms)
  return(out)
}


# Model -------------------------------------------------------------------

modelPasin <- function(t,y,parms){
  with(as.list(c(y, parms)), {
    dA = phi_S * exp(-delta_S*t)+phi_L * exp(-delta_L*t)-delta_Ab*A

    list(dA)
  })
}

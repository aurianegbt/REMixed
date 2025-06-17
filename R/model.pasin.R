#' Model from Pasin and al.,2019
#'
#' @description
#' Generate trajectory of the Humoral Immune Response to a Prime-Boost Ebola Vaccine.
#'
#' @details
#' The model correspond to the dynamics of the humoral response, from 7 days after the boost immunization with antibodies secreting cells -\eqn{S} and \eqn{L}, characterized by their half lives- that produces antibodies -\eqn{AB}- at rate \eqn{\theta_S} and \eqn{\theta_L}. All these biological entities decay at rate repectively \eqn{\delta_S, \delta_L} and \eqn{\delta_{Ab}}. Model is then defined as
#' \deqn{ \left\{\begin{matrix}\frac{d}{dt} Ab(t) &=& \theta_S S(t) + \theta_L L(t) - \delta_{Ab} Ab(t)  \\ \frac{d}{dt}S(t) &=& -\delta_S S(t) \\ \frac{d}{dt} L(t) &=& -\delta_L L(t)\end{matrix}\right. }
#'
#' @param t vector of time ;
#' @param y initial condition, named vector of form \code{c(Ab=<...>,S=<...>,L=<...>)} ;
#' @param parms named vector of model parameter ; should contain "\code{theta_S}","\code{theta_L}","\code{delta_Ab}","\code{delta_S}","\code{delta_L}".
#'
#' @return Matrix of \code{time} and observation of antibody titer \code{Ab}, and ASCs \code{S} and \code{L}.
#' @export
#' @seealso \code{\link{indParm}}, \code{\link{model.clairon}}.
#'
#' @examples
#' y = c(Ab=0,S=5,L=5)
#' parms = c(theta_S = 611,
#'           theta_L = 3.5,
#'           delta_Ab = 0.025,
#'           delta_S = 0.231,
#'           delta_L = 0.000152)
#'
#' t = seq(0,100,5)
#' res <- model.pasin(t,y,parms)
#' plot(res)
#'
#' @references Pasin C, Balelli I, Van Effelterre T, Bockstal V, Solforosi L, Prague M, Douoguih M, Thiébaut R, for the EBOVAC1 Consortium. 2019. Dynamics of the humoral immune response to a prime-boost Ebola vaccine: quantification and sources of variation. J Virol 93: e00579-19. https://doi.org/10.1128/JVI.00579-19


model.pasin <- function(t,y,parms){
  #phi_S,phi_L,delta_Ab,delta_S=0.23,delta_L=0.000316
  if(!setequal(names(y),c("Ab","S","L"))){
    bonus = setdiff(names(y),c("Ab","L","S"))
    malus = setdiff(c("Ab","S","L"),names(y))
    if(length(bonus)!=0 & length(malus) !=0){
    stop(paste0("Missing initial condition for ",malus," and ",bonus," isn't in the model."))
    }else if(length(bonus)!=0){
      stop(paste0(bonus," isn't dynamic of the model."))
    }else if(length(malus) !=0){
      stop(paste0("Missing initial condition for ",malus))
    }
  }
  if(!setequal(names(parms),c("theta_S","theta_L","delta_Ab","delta_S","delta_L"))){
    stop(paste0("Missing parmeters ",setdiff(c("theta_S","theta_L","delta_Ab","delta_S","delta_L"),names(parms))," and ",setdiff(names(parms),c("theta_S","theta_L","delta_Ab","delta_S","delta_L"))," isn't in the model."))
  }
  parms <- c(parms,S0=y[["S"]],L0=y[["L"]])
  out <- deSolve::ode(y["Ab"],t,
                      function(t,y,parms){
                        with(as.list(c(y, parms)), {
                          dAb = theta_S * S0 * exp(-delta_S*t) + theta_L * L0* exp(-delta_L*t)  - delta_Ab * Ab
                          list(c(dAb))
                        })
                      },parms)
  out <- cbind(out,
               S=y[["S"]]*exp(-parms[["delta_S"]]*out[,"time"]),
               L=y[["L"]]*exp(-parms[["delta_L"]]*out[,"time"]))
  return(out)
}


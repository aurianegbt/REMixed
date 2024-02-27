#' Antibodies concentration
#'
#' @param t time
#' @param A0 antibody at baseline
#' @param phi_S phi_S
#' @param phi_L phi_L
#' @param delta_Ab delta_Ab
#' @param delta_S delta_S
#' @param delta_L delta_L
#'
#' @return éa
#' @export
#'
#' @examples a
Pasin <- function(t,A0,phi_S,phi_L,delta_Ab,delta_S=0.23,delta_L=0.000316){

  y =c(A=A0)
  parms=c(phi_S=phi_S,
          phi_L=phi_L,
          delta_Ab=delta_Ab,
          delta_S=delta_S,
          delta_L=delta_L)

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

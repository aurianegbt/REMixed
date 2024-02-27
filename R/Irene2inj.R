#' Title
#'
#' @param t time
#' @param S0 s0
#' @param A0 a
#' @param t_inj a
#' @param fM1 a
#' @param fM2 a
#' @param theta1 a
#' @param theta_2 a
#' @param delta_S a
#' @param delta_V a
#' @param delta_Ab a
#'
#' @return a
#' @export
#'
#' @examples a
Covid <- function(t,S0,A0,t_inj,fM1,fM2,theta1,theta2,delta_S,delta_V=2.7,delta_Ab=0.03){
  y=c(S=S0,A=A0)
  parms=c(t_0=t[1],
          t_inj=t_inj,
          fM1=fM1,
          fM2=fM2,
          theta1=theta1,
          theta2=theta2,
          delta_S=delta_S,
          delta_Ab=delta_Ab,
          delta_V=delta_V)

  out <- deSolve::ode(y,t,modelIrene2inj,parms)
  return(out)
}


# Model  ------------------------------------------------------------------

modelIrene2inj <- function(t,y,parms){
  with(as.list(c(y, parms)), {

    if(t < t_inj){
      C = fM1
      ttilde = t_0
      theta = theta1
    }else{
      C = fM2
      ttilde = t_inj
      theta = theta2
    }

    dS = C*exp(-delta_V*(t-ttilde))-delta_S*S
    dA = theta * S - delta_Ab * A

    list(c(dS,dA))
  })
}

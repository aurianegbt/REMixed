#' Generate trajectory of antibody response to SARS-Cov vaccine.
#'
#' The model correspond to the antibody response to two injection of vaccine at time `t_0`and `t_inj`.Model is defined as
#' \deqn{\displaystyle\left\{\begin{matrix} \frac{d}{dt} S(t) &=& f_{\overline M_k} e^{-\delta_V(t-t_k)}-\delta_S S(t) \\ \frac{d}{dt} Ab(t) &=& \theta S(t) - \delta_{Ab} Ab(t)\end{matrix}\right.  }
#' on each interval \eqn{I_1=[t_0;t_{inj}[ } and \eqn{I_2=[t_{inj};+\infty[}. For each interval \eqn{I_k}, we have \eqn{t_k} corresponding to the last injection date of vaccine, such that \eqn{t_1=t_0} and \eqn{t_2=t_{inj}}. By definition, \eqn{f_{\overline M_1}=1} (Clairon and al., 2023).
#'
#' @param t vector of time ;
#' @param y initial condition, named vector of form c(S=S0,A=A0) ;
#' @param parms named vector of model parameter ; should contain "t_0","t_inj","fM2","theta","delta_S","delta_Ab","delta_V".
#'
#' @return Matrix of time and observation of antibody titer Ab.
#' @export
#'
#' @examples
#' y = c(S=1,A=0)
#'
#' parms = c(t_0 = 0,
#'           t_inj = 21,
#'           fM2 = 4.5,
#'           theta = 18.7,
#'           delta_S = 0.01,
#'           delta_Ab = 0.23,
#'           delta_V = 2.7)
#'
#' t = seq(0,35,1)
#'
#' res <- Clairon(t,y,parms)
#'
#' plot(res)
#'
#' @references Quentin Clairon, Mélanie Prague, Delphine Planas, Timothée Bruel, Laurent Hocqueloux, et al.. Modeling the evolution of the neutralizing antibody response against SARS-CoV-2 variants after several administrations of Bnt162b2. 2023. hal-03946556

Clairon <- function(t,y,parms){
  #t_inj,fM1,fM2,theta1,theta2,delta_S,delta_V=2.7,delta_Ab=0.03
  if(!setequal(names(y),c("S","A"))){
    stop(paste0("Missing initial condition for ",setdiff(c("S","A"),names(y))," and ",setdiff(names(y),c("S","A"))," isn't in the model."))
  }
  if(!setequal(names(parms),c("t_0","t_inj","fM2","theta","delta_S","delta_Ab","delta_V"))){
    stop(paste0("Missing parmeters ",setdiff(c("t_0","t_inj","fM2","theta","delta_S","delta_Ab","delta_V"),names(parms))," and ",setdiff(names(parms),c("t_0","t_inj","fM2","theta","delta_S","delta_Ab","delta_V"))," isn't in the model."))
  }
  out <- deSolve::ode(y,t,modelClairon,parms)
  return(out)
}


# Model  ------------------------------------------------------------------

modelClairon <- function(t,y,parms){
  with(as.list(c(y, parms)), {

    if(t < t_inj){
      C = 1
      ttilde = t_0
    }else{
      C = fM2
      ttilde = t_inj
    }

    dS = C*exp(-delta_V*(t-ttilde))-delta_S*S
    dA = theta * S - delta_Ab * A

    list(c(dS,dA))
  })
}

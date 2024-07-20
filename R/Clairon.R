#' Generate trajectory of antibody response to SARS-Cov vaccine.
#'
#' The model correspond to the antibody response to two injection of vaccine at time \eqn{t_0=0} and \eqn{t_{inj}}.Model is defined as
#' \deqn{\displaystyle\left\{\begin{matrix} \frac{d}{dt} S(t) &=& f_{\overline M_k} e^{-\delta_V(t-t_k)}-\delta_S S(t) \\ \frac{d}{dt} AB(t) &=& \theta S(t) - \delta_{Ab} AB(t)\end{matrix}\right.  }
#' on each interval \eqn{I_1=[0;t_{inj}[ } and \eqn{I_2=[t_{inj};+\infty[}. For each interval \eqn{I_k}, we have \eqn{t_k} corresponding to the last injection date of vaccine, such that \eqn{t_1=0} and \eqn{t_2=t_{inj}}. By definition, \eqn{f_{\overline M_1}=1} (Clairon and al., 2023).
#'
#' @param t vector of time ;
#' @param y initial condition, named vector of form c(S=S0,AB=A0) ;
#' @param parms named vector of model parameter ; should contain "fM2","theta","delta_S","delta_Ab","delta_V".
#'
#' @return Matrix of time and observation of antibody secreting cells S and antibody titer AB.
#' @export
#' @seealso \code{\link{indParm}}
#'
#' @examples
#' y = c(S=1,AB=0)
#'
#' parms = c(fM2 = 4.5,
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

Clairon <- function(t,y,parms,tinj=21){
  #t_inj,fM1,fM2,theta1,theta2,delta_S,delta_V=2.7,delta_Ab=0.03
  if(!setequal(names(y),c("S","Ab"))){
    stop(paste0("Missing initial condition for ",setdiff(c("S","Ab"),names(y))," and ",setdiff(names(y),c("S","Ab"))," isn't in the model."))
  }
  if(!setequal(names(parms),c("fM2","theta","delta_S","delta_Ab","delta_V"))){
    stop(paste0("Missing parmeters ",setdiff(c("fM2","theta","delta_S","delta_Ab","delta_V"),names(parms))," and ",setdiff(names(parms),c("fM2","theta","delta_S","delta_Ab","delta_V"))," isn't in the model."))
  }
  out <- deSolve::ode(y,t,
                      function(t,y,parms,t_inj=tinj){
                        with(as.list(c(y, parms)), {

                          t_0 = 0

                          if(t < t_inj){
                            C = 1
                            ttilde = t_0
                          }else{
                            C = fM2
                            ttilde = t_inj
                          }

                          dS = C*exp(-delta_V*(t-ttilde))-delta_S*S
                          dAb = theta * S - delta_Ab * Ab

                          list(c(dS,dAb))
                        })
                      },parms)
  return(out)
}

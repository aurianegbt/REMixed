#' Generate trajectory of PK model
#'
#' The administration is via a bolus. The PK model has one compartment (volume V) and a linear elimination (clearance Cl). The parameter ka is defined as \eqn{ ka=\frac{Cl}{V}}.
#'
#' @param t vector of time ;
#' @param y initial condition, named vector of form c(C=C0) ;
#' @param parms named vector of model parameter ; should contain either "Cl" and "V" or "ka".
#'
#' @return Matrix of time and observation of Concentration C.
#' @export
#'
#' @examples
#' res <- PK(seq(0,30,1),c(C=100),parms=c(ka=1))
#'
#' plot(res)

PK <- function(t,y,parms){
  if(!setequal(names(y),"C")){
    stop(paste0("Missing initial condition for ",setdiff(c("C"),names(y))," and ",setdiff(names(y),c("C"))," isn't in the model."))
  }
  if(!setequal(names(parms),c("Cl","V")) & !setequal(names(parms),c("ka")) & setequal(names(parms),c("Cl","V","ka"))){
    stop(paste0("Missing parmeters, parms should include either Cl and V, or ka."))
  }

  if(setequal(names(parms),c("Cl","V"))){
    ka = parms[["Cl"]]/parms[["V"]]
  }else if(setequal(names(parms),c("ka"))){
    ka = parms[["ka"]]
  }else{
    Cl = parms[["Cl"]]
    V = parms[["V"]]
    ka = parms[["ka"]]
    if(ka!=Cl/V){
      warning(paste0("Contradictory values given, Cl and V parameters are ignored and ka=",ka))
    }
  }

  C0 = y[["C"]]

  return(as.matrix(data.frame(time=t,C=C0*exp(-ka*t))))
}

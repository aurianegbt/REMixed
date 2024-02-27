#' Individual PK trajectory
#'
#' @description
#' The administration is via a bolus.
#' The PK model has one compartment (volume $V$) and a linear elimination (clearance $Cl$).
#'
#' @param t vector of observation times.
#' @param Cl Cleareance of subject.
#' @param V Compartment volume of subject.
#' @param C0 Initial concentration $C(t=0)$.
#' @param ka ka value
#'
#' @return The concentration $C(t)$ of subject defined by individual parameters $Cl$, $V$ and $C_0$.
#' @export
#'
#' @examples a NE VA PAS NIVEAU PLOT TIME !!!! TIME A RETOURNER ! VOIR OBJET DESOLVE

PK <- function(t,C0,Cl=NULL,V=NULL,ka=Cl/V){
  if(identical(ka,numeric(0))){
    stop("ka, or, Cl and V, need to be provided.")
  }else if(!is.null(Cl) && !is.null(V) && ka!=Cl/V){
    warning(paste0("Contradictory values given, Cl and V parameters are ignored and ka=",ka))
  }
  return(C0 * exp(-ka*t))
}

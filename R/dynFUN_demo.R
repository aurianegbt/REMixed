#' Dynamic functions demo
#'
#' Example of solver for \code{\link{remix}} and \code{\link{cv.remix}} algorithm. It is perfectly adapted for the Monolix demo project (see \code{\link{getMLXdir}}).
#'
#' @details
#' Suppose you have antibodies secreting cells -\eqn{S}- that produces antibodies -\eqn{AB}- at rate \eqn{\varphi_S}. These two biological entities decay respectively at rate \eqn{\delta_S} and \eqn{\delta_{AB}}. The biological mechanism behind is :
#' \deqn{ \left\{\begin{matrix} \frac{d}{dt}S(t) &=& -\delta_S S(t) \\ \frac{d}{dt} AB(t) &=& \varphi_S S(t) - \delta_{AB} AB(t)  \\ (S(0),AB(0)) &=& (S_0,AB_0) \end{matrix}\right. }
#'
#' @format \code{dynFUN_demo}
#' function of \code{t}, \code{y}, \code{parms} :
#' \describe{
#'   \item{\code{t}}{vector of timepoint.}
#'   \item{\code{y}}{initial condition, named vector of form \code{c(AB=<...>,S=<...>)}.}
#'   \item{\code{parms}}{named vector of model parameter ; should contain \code{phi_S},\code{delta_AB},\code{delta_S}.}
#' }
#'
#' @seealso \code{\link{model.pasin}}, \code{\link{getMLXdir}}.
#'
#' @examples
#' t = seq(0,300,1)
#' y =c(AB=1000,S=5)
#' parms = c(phi_S = 611, delta_AB = 0.03, delta_S=0.01)
#'
#' res <- dynFUN_demo(t,y,parms)
#'
#' plot(res[,"time"],
#'      log10(res[,"AB"]),
#'      ylab="log10(AB(t))",
#'      xlab="time (days)",
#'      main="Antibody titer over the time",
#'      type="l")
#'
#' plot(res[,"time"],
#'      res[,"S"],
#'      ylab="S(t)",
#'      xlab="time (days)",
#'      main="Antibody secreting cells quantity over time",
#'      type="l")
#'
#' @references Pasin C, Balelli I, Van Effelterre T, Bockstal V, Solforosi L, Prague M, Douoguih M, Thiébaut R, for the EBOVAC1 Consortium. 2019. Dynamics of the humoral immune response to a prime-boost Ebola vaccine: quantification and sources of variation. J Virol 93 : e00579-19. https://doi.org/10.1128/JVI.00579-19

"dynFUN_demo"

#' Remix results
#'
#' Extracts the build minimizing the BIC over a grid of lambda.
#'
#'
#' @param fit output of \code{\link{cv.Remix}}.
#'
#' @return outputs from \code{\link{Remix}} algorithm achieving the best BIC among those computed by \code{\link{cv.Remix}}.
#' @export
#' @seealso \code{\link{cv.Remix}}, \code{\link{Remix}}.
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(AB=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' y = c(S=5,AB=1000)
#'
#' cv.outputs = cv.Remix(project = project,
#'             dynFUN = dynFUN_demo,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             ncores=8,
#'             eps2=1)
#'
#' res <- retrieveBest(cv.outputs)
#'
#' plotConvergence(res)
#'
#' trueValue = read.csv(paste0(dirname(project),"/demoSMLX/Simulation/populationParameters.txt"))#'
#'
#' plotSAEM(res,paramToPlot = c("delta_S_pop","phi_S_pop","delta_AB_pop"),trueValue=trueValue)
#' }
retrieveBest <- function(fit){
  if(!inherits(fit,"cvRemix")){
    stop("")
  }

  argmin.id = which.min(fit$BIC)

  results = list(info = append(fit$info,lambda=fit$lambda[argmin.id]),
                 finalRes = fit$res[[argmin.id]],
                 iterOutputs = fit$outputs[[argmin.id]])
  class(results) <- "remix"
  return(results)
}

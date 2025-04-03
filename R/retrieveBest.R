#' Remix results
#'
#' Extracts the build minimizing an information criterion over a grid of lambda.
#'
#'
#' @param fit output of \code{\link{cv.remix}};
#' @param criterion which criterion function to take into account. Default is the function 'BIC", but one can use 'eBIC', 'AIC', 'AICc' or any function depending on a `cvRemix` object.
#'
#' @return outputs from \code{\link{remix}} algorithm achieving the best IC among those computed by \code{\link{cv.remix}}.
#' @export
#' @seealso \code{\link{cv.remix}}, \code{\link{remix}}, \code{\link{BIC.cvRemix}}, \code{\link{eBIC.cvRemix}}, \code{\link{AIC.cvRemix}}, \code{\link{AICc.cvRemix}}.
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
retrieveBest <- function(fit,criterion=BIC){
  if(!inherits(fit,"cvRemix")){
    stop("Class of fit must be cvRemix")
  }
  argmin.id = which.min(criterion(fit))

  results = list(info = append(fit$info,list(lambda=fit$lambda[argmin.id],finalSAEM=FALSE,test=FALSE,p.max=NULL,project=if(length(fit$info$project)==1){fit$info$project}else{fit$info$project[argmin.id]}))[c("project","param.toprint","regParam.toprint","alpha","lambda","finalSAEM","test","p.max","N","ntot")],
                 finalRes = append(fit$res[[argmin.id]],list(saemBeforeTest=NULL)),
                 iterOutputs = fit$outputs[[argmin.id]])
  class(results) <- "remix"
  return(results)
}

#' extract remix results from cvRemix object
#'
#' Extracts a build from a cvRemix object.
#'
#'
#' @param fit output of \code{\link{cv.remix}};
#' @param n rank (in the `fit$lambda`) to extract.
#'
#' @return outputs from \code{\link{remix}} algorithm of rank `n` computed by \code{\link{cv.remix}}.
#' @export
#' @seealso \code{\link{cv.remix}}, \code{\link{remix}}, \code{\link{BIC.cvRemix}}, \code{\link{eBIC.cvRemix}}, \code{\link{AIC.cvRemix}}, \code{\link{AICc.cvRemix}}.
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
#' res <- extract(cv.outputs,6)
#'
#' plotConvergence(res)
#'
#' trueValue = read.csv(paste0(dirname(project),"/demoSMLX/Simulation/populationParameters.txt"))
#'
#'
#' plotSAEM(res,paramToPlot = c("delta_S_pop","phi_S_pop","delta_AB_pop"),trueValue=trueValue)
#' }
extract <- function(fit,n){
  if(!inherits(fit,"cvRemix")){
    stop("Class of fit must be cvRemix")
  }
  if(!(n>=1 && n<=length(fit$lambda))){
    stop(paste0("n must be between 1 and ",lenth(fit$lambda)))
  }

  results = list(info = append(fit$info,list(lambda=fit$lambda[n],finalSAEM=FALSE,test=FALSE,p.max=NULL,project=if(length(fit$info$project)==1){fit$info$project}else{fit$info$project[n]}))[c("project","param.toprint","regParam.toprint","alpha","lambda","finalSAEM","test","p.max","N","ntot")],
                 finalRes = append(fit$res[[n]],list(saemBeforeTest=NULL)),
                 iterOutputs = fit$outputs[[n]])
  class(results) <- "remix"
  return(results)
}

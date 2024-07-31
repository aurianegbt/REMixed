#' Display the value of parameters at each iteration
#'
#' @param fit object of class remix, from \code{\link{Remix}} or a certain build from \code{\link{cv.Remix}} output.
#' @param paramToPlot Population parameters to plot (which have been estimated by SAEM) ;
#' @param trueValue (for simulation purpose) vector named of true values ;
#'
#' @return
#' For each parameters, the values at the end of each iteration of remix algorithm is drawn. Moreover, the SAEM steps of each iteration are displayed.
#' @export
#'
#' @seealso \code{\link{Remix}}, \code{\link{cv.Remix}}.
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
#' lambda = 1440
#'
#' res = Remix(project = project,
#'             dynFUN = dynFUN_demo,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             eps2=1,
#'             lambda=lambda)
#'
#' plotConvergence(res)
#'
#' trueValue = read.csv(paste0(dirname(project),"/demoSMLX/Simulation/populationParameters.txt"))#'
#'
#' plotSAEM(res,paramToPlot = c("delta_S_pop","phi_S_pop","delta_AB_pop"),trueValue=trueValue)
#' }
plotSAEM <- function(fit,paramToPlot = 'all',trueValue=NULL){ # GENES are GENES to remove (id with name of variable )
  estimates = dplyr::filter(fit$iterOutputs$estimates,iteration!=0)
  if(!identical(paramToPlot,"all")){
    estimates = dplyr::select(estimates,iteration, phase, all_of(paramToPlot))
  }

  ## Repérer les changements d'itération, i.e. phase de 2 à 1 # "#e0e6c6","#cac2ba","#ffe591"
  phases <- data.frame(start=1,end=1,phase="Initialisation",estimates[1,paramToPlot])
  start = 1
  prev = 1
  remix.iteration = 0
  for(i in 1:(nrow(estimates)-1)){
    new = estimates$phase[i+1]
    if(prev==2 && new==1){

      phases <- rbind(phases,
                      cbind(data.frame(start = start, end= i, phase = ifelse(remix.iteration==0,"Initialisation",paste0("Iteration ",remix.iteration))), estimates[i,paramToPlot]))

      remix.iteration <- remix.iteration+1
      start <- i
    }
    prev <- new
  }
  phases <- rbind(phases,
                  cbind(data.frame(start = start, end= nrow(estimates), phase = ifelse(remix.iteration==0,"Initialisation",paste0("Iteration ",remix.iteration))),estimates[nrow(estimates),paramToPlot]))
  color <- c("#cac2ba",rep(c("#e0e6c6","#ffe591"),nrow(phases)))[1:length(unique(phases$phase))]

  estimates$iteration <- 1:nrow(estimates)

  dim = ceiling(sqrt(length(paramToPlot)))

  for(p in paramToPlot){
    cmd <- paste0(" plot.",p," <- ggplot2::ggplot() +
    ggplot2::geom_line(data=estimates,ggplot2::aes(x=iteration,y=",p,"),color='grey38') +
    ggplot2::geom_rect(data=phases,ggplot2::aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf,fill=phase),alpha=0.3)+
    ggplot2::geom_line(data=phases,ggplot2::aes(x=end,y=",p,"),color='indianred',lwd=1)+
    ggplot2::geom_point(data=phases,ggplot2::aes(x=end,y=",p,"),color='indianred',size=2)+
    ggplot2::scale_fill_manual(values=color,breaks=phases$phase[-1]) +
    ggplot2::theme(legend.position='none')+
    ggplot2::xlim(c(0,nrow(estimates)))")
    if(!is.null(trueValue)){
      cmd <- paste0(cmd,paste0("+
    ggplot2::geom_segment(ggplot2::aes(x=0,xend=+Inf,y=trueValue[['",p,"']],yend=trueValue[['",p,"']]),color='slateblue',lwd=1)"))
    }
    eval(parse(text=cmd))
  }

  cmd <- paste0("finalPlot <-     ggpubr::annotate_figure(ggpubr::ggarrange(",paste0('plot.',paramToPlot,collapse=", "),", ncol=dim,nrow=dim),top=ggpubr::text_grob('Values of parameter throughout each SAEM iterations',face='bold',size=18,color='#862B0D'))")

  eval(parse(text=cmd))

  return(finalPlot)
}

#' Log-likelihood convergence.
#'
#' @param fit fit object of class remix, from \code{\link{Remix}} or a certain build from \code{\link{cv.Remix}} output.
#'
#' @return Log-Likelihood values throughout the algorithm iteration.
#' @export
#'
#' @seealso \code{\link{Remix}}, \code{\link{cv.Remix}}.
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
#' lambda = 1440
#'
#' res = Remix(project = project,
#'             dynFUN = dynFUN,
#'             y = y,
#'             ObsModel.transfo = ObsModel.transfo,
#'             alpha = alpha,
#'             selfInit = TRUE,
#'             eps1=10**(-2),
#'             eps2=1,
#'             lambda=lambda)
#'
#' plotConvergence(res)
#'
#' trueValue = read.csv(paste0(dirname(project),"/demoSMLX/Simulation/populationParameters.txt"))#'
#'
#' plotSAEM(res,paramToPlot = c("delta_S_pop","phi_S_pop","delta_AB_pop"),trueValue=trueValue)
#' }
plotConvergence <- function(fit){
  LL.outputs <- sapply(fit$iterOutputs$LL,FUN=function(Liter){Liter$LL})
  LLpen.outputs <- unlist(fit$iterOutputs$LL.pen)

  if(length(LL.outputs) != length(LLpen.outputs)){
    stop("Different number of outputs for log-likelihood and penalised log-likelihood.")
  }

  dfToPlot <- data.frame(iteration = 1:length(LL.outputs),LL=LL.outputs,LLpen=LLpen.outputs)


  plotUp <-
    ggplot2::ggplot(dfToPlot,ggplot2::aes(x=iteration,y=LLpen))+
      ggplot2::geom_line()+
      ggplot2::geom_point() +
      ggplot2::ggtitle("Penalised log-likelihood value throughout the iteration. ") +
      ggplot2::ylab("Penalised log-likelihood")

  return(invisible(plotUp))
}

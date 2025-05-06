#' Plot of cv.remix object
#'
#' Calibration and Informative Criterion plot for cvRemix object.
#'
#' @param x output of \code{\link{cv.remix}}.
#'
#' @return A plot.
#'
#' @seealso \code{\link{cv.remix}}
#'
#' @export
plot.cvRemix <- function(x,criterion=BICc,trueValue=NULL){
  pC <- plotCalibration(x,criterion=criterion,trueValue = trueValue)+ggplot2::ggtitle("")
  pI <- plotIC(x,criterion=criterion)+ggplot2::ggtitle("")

  plot <-
    ggpubr::ggarrange(pC,pI,
                      labels = c("Calibration Plot          ",
                                 "Information Criterion Plot"),
                      hjust=-0.2,
                      font.label = list(size = 14, color = "brown4", face = "bold"),
                      ncol=1)

  return(plot)
}

#' Display the value of parameters at each iteration
#'
#' @param fit object of class remix, from \code{\link{remix}} or a certain build from \code{\link{cv.remix}} output.
#' @param paramToPlot Population parameters to plot (which have been estimated by SAEM) ;
#' @param trueValue (for simulation purpose) vector named of true values ;
#'
#' @return
#' For each parameters, the values at the end of each iteration of remix algorithm is drawn. Moreover, the SAEM steps of each iteration are displayed.
#' @export
#'
#' @seealso \code{\link{remix}}, \code{\link{cv.remix}}.
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
#' res = remix(project = project,
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
#' trueValue = read.csv(paste0(dirname(project),"/demoSMLX/Simulation/populationParameters.txt"))
#'
#' plotSAEM(res,paramToPlot = c("delta_S_pop","phi_S_pop","delta_AB_pop"),trueValue=trueValue)
#' }
plotSAEM <- function(fit,paramToPlot = 'all',trueValue=NULL){ # GENES are GENES to remove (id with name of variable )
  if(!inherits(fit,"remix")){
    stop("Class of fit must be remix")
  }
  if(identical(paramToPlot,"all")){
    paramToPlot <- fit$info$param.toprint
  }
  iteration <- phase  <- finalPlot <-  NULL
  estimates = dplyr::filter(fit$iterOutputs$estimates,iteration!=0)
  if(!identical(paramToPlot,"all")){
    estimates = dplyr::select(estimates,iteration, phase, dplyr::all_of(paramToPlot))
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

  finalPlot
}

#' Log-likelihood convergence.
#'
#' @param fit fit object of class remix, from \code{\link{remix}} or a certain build from \code{\link{cv.remix}} output.
#'
#' @return Log-Likelihood values throughout the algorithm iteration.
#' @export
#'
#' @seealso \code{\link{remix}}, \code{\link{cv.remix}}.
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
#' res = remix(project = project,
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
  if(!inherits(fit,"remix")){
    stop("Class of fit must be remix")
  }


  iteration <-  NULL

  niter = fit$finalRes$iter
  LL.outputs <- sapply(fit$iterOutputs$LL,FUN=function(Liter){Liter$LL})
  LL.pen <- sapply(1:niter,FUN=function(i){
    LL.outputs[i] - fit$info$lambda*sum(fit$iterOutputs$param[i,fit$info$regParam.toprint]!=0)
  })

  dfToPlot <- data.frame(iteration = 1:length(LL.outputs),LL=LL.outputs,LL.pen = LL.pen)


  ggplot2::ggplot(dfToPlot,ggplot2::aes(x=iteration,y=LL.pen))+
    ggplot2::geom_line()+
    ggplot2::geom_point() +
    ggplot2::ggtitle("Penalised log-likelihood value throughout the iteration. ") +
    ggplot2::ylab("Penalised log-likelihood")
}



#' Calibration plot.
#'
#' @param fit fit object of class cvRemix, from \code{\link{cv.remix}}.
#' @param legend.position (default NULL) 	the default position of legends ("none", "left", "right", "bottom", "top", "inside").
#' @param trueValue (for simulation purpose) named vector containing the true value of regularization parameter.
#' @param criterion function ; which criterion among 'BIC', 'eBIC', 'AIC', 'BICc', or function of cvRemix object to take into account (default : BICc).
#' @param dismin logical ; if minimizers of information criterion should be display.
#'
#' @return Calibration plot, over the lambda.grid.
#' @export
#'
#' @seealso \code{\link{remix}}, \code{\link{cv.remix}}.
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
#' res = cv.remix(project = project,
#'                dynFUN = dynFUN_demo,
#'                y = y,
#'                ObsModel.transfo = ObsModel.transfo,
#'                alpha = alpha,
#'                selfInit = TRUE,
#'                eps1=10**(-2),
#'                ncores=8,
#'                nlambda=8,
#'                eps2=1)
#'
#' plotCalibration(res)
#'
#' plotIC(res)
#' }
plotCalibration <- function(fit,legend.position = "none",trueValue=NULL,criterion=BICc,dismin=TRUE){
  if(!inherits(fit,"cvRemix")){
    stop("Class of fit must be cvRemix")
  }

  null <- rep(0,length(fit$info$alpha$alpha1))
  null[as.logical(trueValue[fit$info$alpha$alpha1]!=0)] <- 1:sum(trueValue[fit$info$alpha$alpha1]!=0)

  df <- data.frame()

  for(i in 1:length(fit$lambda)){
    df <- rbind(df,data.frame(parameter = unname(fit$info$alpha$alpha1),
                              value = fit$res[[i]]$alpha,
                              lambda = fit$lambda[[i]],
                              criterion = criterion(extract(fit,i)),
                              LL = fit$res[[i]]$LL))

  }

  if(!is.null(trueValue)){
    trueValueDF <- data.frame(parameter=fit$info$alpha$alpha1[as.logical(trueValue[fit$info$alpha$alpha1]!=0)],
                              y=as.numeric(trueValue[fit$info$alpha$alpha1[as.logical(trueValue[fit$info$alpha$alpha1]!=0)]]),
                              yend=as.numeric(trueValue[fit$info$alpha$alpha1[as.logical(trueValue[fit$info$alpha$alpha1]!=0)]]),
                              null=1:sum(trueValue[fit$info$alpha$alpha1]!=0))

    df <- cbind(df,null = rep(null,length(fit$lambda)))

    p <-
      ggplot2::ggplot(df,ggplot2::aes(x=lambda,y=value,color=as.factor(null),group=parameter)) +
      ggplot2::theme(legend.position = legend.position) +
      ggplot2::geom_line() +
      ggplot2::geom_point()  +
    ggplot2::geom_segment(data=trueValueDF,mapping=ggplot2::aes(color=factor(null),group=parameter,y=y,yend=yend),x=-Inf,xend=+Inf) +
      ggplot2::xlab("\u03bb")+
      ggplot2::ylab("\u03b1")
    if(dismin){
      p <- p +
        ggplot2::geom_segment(x = df[which.min(df$criterion),"lambda"],xend = df[which.min(df$criterion),"lambda"], y = -Inf, yend=+Inf,col="indianred",lwd=1)
    }
  }else{
    p <-
      ggplot2::ggplot(df,ggplot2::aes(x=lambda,y=value,color=parameter,group=parameter)) +
      ggplot2::theme(legend.position = legend.position) +
      ggplot2::geom_line() +
      ggplot2::geom_point()  +
      ggplot2::xlab("\u03bb")+
      ggplot2::ylab("\u03b1")
    if(dismin){
      p <- p +
        ggplot2::geom_segment(x = df[which.min(df$criterion),"lambda"],xend = df[which.min(df$criterion),"lambda"], y = -Inf, yend=+Inf,col="indianred",lwd=1)
    }
  }
  return(p)
}

#' IC plot.
#'
#' @param fit fit object of class cvRemix, from \code{\link{cv.remix}};
#' @param criterion which criterion among 'BICc', 'BIC', 'AIC' or 'eBIC' to take into account (default: BICc);
#' @param dismin logical ; if minimizers of information criterion should be display.
#'
#' @return IC trhoughout the lambda.grid.
#' @export
#'
#' @seealso \code{\link{remix}}, \code{\link{cv.remix}}.
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
#' res = cv.remix(project = project,
#'                dynFUN = dynFUN_demo,
#'                y = y,
#'                ObsModel.transfo = ObsModel.transfo,
#'                alpha = alpha,
#'                selfInit = TRUE,
#'                eps1=10**(-2),
#'                ncores=8,
#'                nlambda=8,
#'                eps2=1)
#'
#' plotCalibration(res)
#'
#' plotIC(res)
#' }
plotIC <- function(fit,criterion=BICc,dismin=TRUE){
  if(!inherits(fit,"cvRemix")){
  stop("Class of fit must be cvRemix")
  }


  df <- data.frame(criterion = sapply(1:length(fit$lambda),FUN=function(i){criterion(extract(fit,i))}),lambda=fit$lambda)


  p <-
  ggplot2::ggplot(df,ggplot2::aes(x=lambda,y=criterion)) +
    ggplot2::geom_line(color="slateblue") + ggplot2::geom_point(color="slateblue4",size=1,shape=3) +
    ggplot2::ggtitle(paste0('Criterion throughout the grid of \u03bb')) +
    ggplot2::xlab('\u03bb') + ggplot2::theme(plot.title = ggplot2::element_text(size=20,color="red4"))
  if(dismin){
    p <- p +
      ggplot2::geom_segment(x = df[which.min(df$criterion),"lambda"],xend = df[which.min(df$criterion),"lambda"], y = -Inf, yend=+Inf,col="indianred",lwd=0.7,linetype=5) +
      ggplot2::geom_segment(y = df[which.min(df$criterion),"criterion"],yend = df[which.min(df$criterion),"criterion"], x = -Inf, xend=+Inf,col="indianred",lwd=0.7,linetype=5) +
      ggplot2::geom_point(y = df[which.min(df$criterion),"criterion"],x = df[which.min(df$criterion),"lambda"],col="indianred4",size=2)+
      ggplot2::geom_text(label = paste0("\u03bb.min = ",round(df[which.min(df$criterion),"lambda"],digits=2),"\n IC.min = ",round(df[which.min(df$criterion),"criterion"],digits=2)), x= df[which.max(df$lambda),"lambda"],y = df[which.min(df$criterion),"criterion"],hjust=1,vjust=-0.2,color="indianred4",size=4)
  }
  return(p)
}

#' Plot initialization
#'
#' @param init outputs from \code{\link{initStrat}} function.
#' @param alpha named list of named vector "\code{alpha0}", "\code{alpha1}" (all \code{alpha1} are mandatory). The name of \code{alpha$alpha0} and \code{alpha$alpha1} are the observation model names from the monolix project to which they are linked (if the observations models are defined whithout intercept, alpha$alpha0 need to be set to the vector NULL).
#' @param trueValue (for simulation purpose) named vector containing the true value of regularization parameter.
#'
#' @returns log-likelihood value for all groups of genes tested.
#' @export
#'
#' @seealso \code{\link{initStrat}}.
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
#' init <- initStrat(project,alpha,ObsModel.transfo,seed=1710)
#'
#' plotInit(init)
#' }
plotInit <- function(init,alpha=NULL,trueValue=NULL){
  if(!inherits(fit,"init")){
    stop("Class of init must be init")
  }

  results <- data.frame(genes = unname(sapply(init,function(i){paste0("(",paste0(sort(i$genes),collapse=","),")")})),
                        LL = sapply(init,function(i){-1/2*i$LL[["OFV"]]}))

  if(is.null(trueValue)){
    ggplot2::ggplot(results,ggplot2::aes(x=genes,y=LL))+ggplot2::geom_point(size=3)+
      ggplot2::ylab("Log-Likelihood")+
      ggplot2::xlab("Genes pair") +
      ggplot2::theme(axis.title=ggplot2::element_text(size=12),
            axis.text.y = ggplot2::element_text(size=10),
            axis.text.x = ggplot2::element_text(size=8,angle=90,vjust = 0.5),
            strip.text = ggplot2::element_text(size = 11),
            plot.background = ggplot2::element_rect(fill='transparent', color=NA),
            legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
            legend.background = ggplot2::element_rect(fill = "transparent", color = NA))
  }else{
    if(!is.null(alpha)){
      stop("alpha arguments need to be provided")
    }
    genes = data.frame(yGk=names(alpha$alpha1),
                       true = as.numeric(trueValue[alpha$alpha1])!=0)

    results <- cbind(results,
                     P = sapply(init,function(i){sum(genes[genes$yGk %in% i$genes,"true"])}))
  }

  ggplot2::ggplot(results,ggplot2::aes(x=genes,y=LL,color=P))+ggplot2::geom_point(size=3)+
    ggplot2::ylab("Log-Likelihood")+
    ggplot2::xlab("Genes pair") +
    ggplot2::theme(axis.title=ggplot2::element_text(size=12),
                   axis.text.y = ggplot2::element_text(size=10),
                   axis.text.x = ggplot2::element_text(size=8,angle=90,vjust = 0.5),
                   strip.text = ggplot2::element_text(size = 11),
                   plot.background = ggplot2::element_rect(fill='transparent', color=NA),
                   legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
                   legend.background = ggplot2::element_rect(fill = "transparent", color = NA)) +
    ggplot2::scale_color_gradient2(name="Contains\nTrue Positives",low ="darkred",high =  "darkseagreen",midpoint=0.5,mid="#999b6e",breaks=0:max(results$P))
}

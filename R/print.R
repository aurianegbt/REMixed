#' @export
print.remix <- function(x,...){
  extra_args <- list(...)

  # Defining default settings
  if (!"trueValue" %in% names(extra_args)) {
    trueValue <- NULL
  }
  if (!"digits" %in% names(extra_args)) {
    digits <- 2
  }

  to.cat <- "\n      - - - <  ESTIMATED PARAMETERS  > - - -     \n\n"
  cat(to.cat)

  if(!is.null(x$finalRes$standardError)){
    sd.est <- x$finalRes$standardError$stochasticApproximation[,-3]
    sd.est <- rbind(sd.est, data.frame(parameter=setdiff(names(x$finalRes$param),sd.est$parameter),se=NA))

    to.print <- data.frame(EstimatedValue = sapply(x$finalRes$param,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)}))

    sd.est <- merge(data.frame(parameter=names(x$finalRes$param),EstimatedValue=unname(x$finalRes$param)),sd.est,by="parameter")
    sd.est <- cbind(sd.est,CI_95=sapply(1:nrow(sd.est),FUN=function(i){ifelse(!is.na(sd.est$se[i]),paste0("[",format(signif(sd.est$EstimatedValue[i]-1.96*sd.est$se[i],digits=digits),scientific=TRUE),";",format(signif(sd.est$EstimatedValue[i]+1.96*sd.est$se[i],digits=digits),scientific=TRUE),"]")," ")}))
    rownames(sd.est) <- sd.est$parameter
    sd.est <- sd.est[rownames(to.print),-1]

    to.print <- dplyr::mutate(dplyr::mutate(sd.est,EstimatedValue = sapply(sd.est$EstimatedValue,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)})),se=sapply(sd.est$se,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)}))
    to.print$se[to.print$se=="NA"] <- " "

    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = format(signif(as.numeric(trueValue[rownames(to.print)]),digits=digits),scientific=TRUE),
                        RelativeBias = round(as.numeric((x$finalRes$param[rownames(to.print)]-trueValue[rownames(to.print)])/trueValue[rownames(to.print)]),digits=digits))
    }

    paramfinal <- x$finalRes$param
    print(to.print)

    to.cat <- "\n      - - - <  CRITERION  > - - -     \n"
    to.cat <- paste0(to.cat,"        LL : ",round(x$finalRes$LL,digits=digits))
    to.cat <- paste0(to.cat,"\n       BIC :  ",BIC(x))
    to.cat <- paste0(to.cat,"\n      BICc :  ",BICc(x),"\n")
    cat(to.cat)
  }else{
    to.print <- data.frame(EstimatedValue = sapply(x$finalRes$param,FUN=function(p){format(signif(p,digits=digits),scientific=TRUE)}))

    if(!is.null(trueValue)){
      to.print <- cbind(to.print,
                        TrueValue = format(signif(as.numeric(trueValue[rownames(to.print)]),digits=digits),scientific=TRUE),
                        RelativeBias = round(as.numeric((x$finalRes$param[rownames(to.print)]-trueValue[rownames(to.print)])/trueValue[rownames(to.print)]),digits=digits))
    }
    print(to.print)

    to.cat <- "\n      - - - <  CRITERION  > - - -     \n"
    to.cat <- paste0(to.cat,"        LL : ",round(x$finalRes$LL,digits=digits))
    to.cat <- paste0(to.cat,"\n       BIC :  ",BIC(x))
    to.cat <- paste0(to.cat,"\n      BICc :  ",BICc(x),"\n")
    cat(to.cat)
  }
}

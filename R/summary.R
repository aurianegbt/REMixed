#' @export
summary.remix <- function(object,...){
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  cat("\n",plain.line)
  cat("Outputs of remix algorithm from the project :\n\t",object$info$project,"\n")
  cat(paste0("with  :\n\t- ",object$info$N," individuals;",
      "\n\t- \u03bb = ",round(object$info$lambda,digits=2),
      if(object$info$finalSAEM){paste0(";\n\t- final SAEM computed",
      if(object$info$test){";\n\t- final test computed"})},
      ".\n"))
  cat(dashed.line)
  cat("\u2022 Final Estimated Parameters :\n")
  if(is.null(object$finalRes$standardError)){
    print(sapply(object$finalRes$param[-which(names(object$finalRes$param) %in% object$info$regParam.toprint)],FUN=function(x){round(x,digits=2)}))
  }else{
    re.param = object$finalRes$param[-which(names(object$finalRes$param) %in% object$info$regParam.toprint)]

    sd.est = object$finalRes$standardError$stochasticApproximation[object$finalRes$standardError$stochasticApproximation$parameter %in% names(re.param),-3]

    to.print <- data.frame(EstimatedValue = sapply(re.param,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    sd.est <- merge(data.frame(parameter=names(re.param),EstimatedValue=unname(re.param)),sd.est,by="parameter")
    sd.est <- cbind(sd.est,CI_95=sapply(1:nrow(sd.est),FUN=function(i){ifelse(!is.na(sd.est$se[i]),paste0("[",format(signif(sd.est$EstimatedValue[i]-1.96*sd.est$se[i],digits=2),scientific=TRUE),";",format(signif(sd.est$EstimatedValue[i]+1.96*sd.est$se[i],digits=2),scientific=TRUE),"]")," ")}))
    rownames(sd.est) <- sd.est$parameter
    sd.est <- sd.est[rownames(to.print),-1]

    to.print <- dplyr::mutate(dplyr::mutate(sd.est,EstimatedValue = sapply(sd.est$EstimatedValue,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)})),se=sapply(sd.est$se,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    print(to.print)
    }
  cat("\t",dashed.short,"\u2022 Final Selected Biomarkers and Estimate:\n")
  if(is.null(object$finalRes$standardError)){
    print(sapply(object$finalRes$param[names(which(object$finalRes$param[object$info$regParam.toprint]!=0))],FUN=function(x){signif(x,digits=3)}))
  }else{

    re.param = object$finalRes$param[names(which(object$finalRes$param[object$info$regParam.toprint]!=0))]

    sd.est = object$finalRes$standardError$stochasticApproximation[object$finalRes$standardError$stochasticApproximation$parameter %in% names(re.param),-3]

    to.print <- data.frame(EstimatedValue = sapply(re.param,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    sd.est <- merge(data.frame(parameter=names(re.param),EstimatedValue=unname(re.param)),sd.est,by="parameter")
    sd.est <- cbind(sd.est,CI_95=sapply(1:nrow(sd.est),FUN=function(i){ifelse(!is.na(sd.est$se[i]),paste0("[",format(signif(sd.est$EstimatedValue[i]-1.96*sd.est$se[i],digits=2),scientific=TRUE),";",format(signif(sd.est$EstimatedValue[i]+1.96*sd.est$se[i],digits=2),scientific=TRUE),"]")," ")}))
    rownames(sd.est) <- sd.est$parameter
    sd.est <- sd.est[rownames(to.print),-1]

    to.print <- dplyr::mutate(dplyr::mutate(sd.est,EstimatedValue = sapply(sd.est$EstimatedValue,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)})),se=sapply(sd.est$se,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    print(to.print)

  }
  cat("\t",dashed.short,"\u2022 Final Scores :\n")
  cat("\t\u00B7 OFV  :",-2*object$finalRes$LL,"\n")
  cat("\t\u00B7 BIC  :",BIC(object),"\n")
  cat("\t\u00B7 BICc :",BICc(object),"\n")
  cat("\nExecuted in",object$finalRes$iter,"iterations, ",round(object$finalRes$time,digits=1),"s.\n")
  cat(plain.line)
}

#' @export
summary.cvRemix <- function(object,...){
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  cat("\n",plain.line)
  cat("Outputs of cv.remix algorithm from the project :\n\t",object$info$project,"\n")
  cat(dashed.line)
  cat("Results over the grid of \u03bb:\n")
  toprint <- apply(data.frame(
                          OFV = -2*object$LL,
                          BIC = BIC(object),
                          BICc = BICc(object),
                          "Positives"=sapply(object$res,FUN=function(r){sum(r$alpha!=0)}),row.names = round(object$lambda,digits=2)),2,FUN=function(x){round(x,digits=2)})
  print(toprint)
  cat(dashed.line)
}

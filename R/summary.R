#' @export
summary.remix <- function(remix.outputs){
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  cat("\n",plain.line)
  cat("Outputs of remix algorithm from the project :\n\t",remix.outputs$info$project,"\n")
  cat("with  :\n\t-",remix.outputs$info$N,"individuals;",
      "\n\t- \u03bb =",round(remix.outputs$info$lambda,digits=2),
      if(remix.outputs$info$finalSAEM){paste0("\n\t- final SAEM computed",
      if(remix.outputs$info$test){";\n\t- final test computed"})},
      ".\n")
  cat(dashed.line)
  cat("\u2022 Final Estimated Parameters :\n")
  if(is.null(remix.outputs$finalRes$standardError)){
    print(sapply(remix.outputs$finalRes$param[-which(names(remix.outputs$finalRes$param) %in% remix.outputs$info$regParam.toprint)],FUN=function(x){round(x,digits=2)}))
  }else{
    re.param = remix.outputs$finalRes$param[-which(names(remix.outputs$finalRes$param) %in% remix.outputs$info$regParam.toprint)]

    sd.est = remix.outputs$finalRes$standardError[remix.outputs$finalRes$standardError$parameter %in% names(re.param),-3]

    to.print <- data.frame(EstimatedValue = sapply(re.param,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    sd.est <- merge(data.frame(parameter=names(re.param),EstimatedValue=unname(re.param)),sd.est,by="parameter")
    sd.est <- cbind(sd.est,CI_95=sapply(1:nrow(sd.est),FUN=function(i){ifelse(!is.na(sd.est$se[i]),paste0("[",format(signif(sd.est$EstimatedValue[i]-1.96*sd.est$se[i],digits=2),scientific=TRUE),";",format(signif(sd.est$EstimatedValue[i]+1.96*sd.est$se[i],digits=2),scientific=TRUE),"]")," ")}))
    rownames(sd.est) <- sd.est$parameter
    sd.est <- sd.est[rownames(to.print),-1]

    to.print <- dplyr::mutate(dplyr::mutate(sd.est,EstimatedValue = sapply(sd.est$EstimatedValue,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)})),se=sapply(sd.est$se,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    print(to.print)
    }
  cat("\t",dashed.short,"\u2022 Final Selected Biomarkers and Estimate:\n")
  if(is.null(remix.outputs$finalRes$standardError)){
    print(sapply(remix.outputs$finalRes$param[names(which(remix.outputs$finalRes$param[remix.outputs$info$regParam.toprint]!=0))],FUN=function(x){signif(x,digits=3)}))
  }else{

    re.param = remix.outputs$finalRes$param[names(which(remix.outputs$finalRes$param[remix.outputs$info$regParam.toprint]!=0))]

    sd.est = remix.outputs$finalRes$standardError[remix.outputs$finalRes$standardError$parameter %in% names(re.param),-3]

    to.print <- data.frame(EstimatedValue = sapply(re.param,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    sd.est <- merge(data.frame(parameter=names(re.param),EstimatedValue=unname(re.param)),sd.est,by="parameter")
    sd.est <- cbind(sd.est,CI_95=sapply(1:nrow(sd.est),FUN=function(i){ifelse(!is.na(sd.est$se[i]),paste0("[",format(signif(sd.est$EstimatedValue[i]-1.96*sd.est$se[i],digits=2),scientific=TRUE),";",format(signif(sd.est$EstimatedValue[i]+1.96*sd.est$se[i],digits=2),scientific=TRUE),"]")," ")}))
    rownames(sd.est) <- sd.est$parameter
    sd.est <- sd.est[rownames(to.print),-1]

    to.print <- dplyr::mutate(dplyr::mutate(sd.est,EstimatedValue = sapply(sd.est$EstimatedValue,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)})),se=sapply(sd.est$se,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    print(to.print)

  }
  cat("\t",dashed.short,"\u2022 Final Scores :\n")
  cat("\t\u00B7 OFV  :",-2*remix.outputs$finalRes$LL,"\n")
  cat("\t\u00B7 BIC  :",BIC(remix.outputs),"\n")
  cat("\t\u00B7 eBIC :",eBIC(remix.outputs),"\n")
  cat("\nExecuted in",remix.outputs$finalRes$iter,"iterations, ",remix.outputs$finalRes$time,"s.\n")
  cat(plain.line)
}

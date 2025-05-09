#' @export
summary.remix <- function(remix.output,...){
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  cat("\n",plain.line)
  cat("Outputs of remix algorithm from the project :\n\t",remix.output$info$project,"\n")
  cat(paste0("with  :\n\t- ",remix.output$info$N," individuals;",
      "\n\t- \u03bb = ",round(remix.output$info$lambda,digits=2),
      if(remix.output$info$finalSAEM){paste0(";\n\t- final SAEM computed",
      if(remix.output$info$test){";\n\t- final test computed"})},
      ".\n"))
  cat(dashed.line)
  cat("\u2022 Final Estimated Parameters :\n")
  if(is.null(remix.output$finalRes$standardError)){
    print(sapply(remix.output$finalRes$param[-which(names(remix.output$finalRes$param) %in% remix.output$info$regParam.toprint)],FUN=function(x){round(x,digits=2)}))
  }else{
    re.param = remix.output$finalRes$param[-which(names(remix.output$finalRes$param) %in% remix.output$info$regParam.toprint)]

    sd.est = remix.output$finalRes$standardError$stochasticApproximation[remix.output$finalRes$standardError$stochasticApproximation$parameter %in% names(re.param),-3]

    to.print <- data.frame(EstimatedValue = sapply(re.param,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    sd.est <- merge(data.frame(parameter=names(re.param),EstimatedValue=unname(re.param)),sd.est,by="parameter")
    sd.est <- cbind(sd.est,CI_95=sapply(1:nrow(sd.est),FUN=function(i){ifelse(!is.na(sd.est$se[i]),paste0("[",format(signif(sd.est$EstimatedValue[i]-1.96*sd.est$se[i],digits=2),scientific=TRUE),";",format(signif(sd.est$EstimatedValue[i]+1.96*sd.est$se[i],digits=2),scientific=TRUE),"]")," ")}))
    rownames(sd.est) <- sd.est$parameter
    sd.est <- sd.est[rownames(to.print),-1]

    to.print <- dplyr::mutate(dplyr::mutate(sd.est,EstimatedValue = sapply(sd.est$EstimatedValue,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)})),se=sapply(sd.est$se,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    print(to.print)
    }
  cat("\t",dashed.short,"\u2022 Final Selected Biomarkers and Estimate:\n")
  if(is.null(remix.output$finalRes$standardError)){
    print(sapply(remix.output$finalRes$param[names(which(remix.output$finalRes$param[remix.output$info$regParam.toprint]!=0))],FUN=function(x){signif(x,digits=3)}))
  }else{

    re.param = remix.output$finalRes$param[names(which(remix.output$finalRes$param[remix.output$info$regParam.toprint]!=0))]

    sd.est = remix.output$finalRes$standardError$stochasticApproximation[remix.output$finalRes$standardError$stochasticApproximation$parameter %in% names(re.param),-3]

    to.print <- data.frame(EstimatedValue = sapply(re.param,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    sd.est <- merge(data.frame(parameter=names(re.param),EstimatedValue=unname(re.param)),sd.est,by="parameter")
    sd.est <- cbind(sd.est,CI_95=sapply(1:nrow(sd.est),FUN=function(i){ifelse(!is.na(sd.est$se[i]),paste0("[",format(signif(sd.est$EstimatedValue[i]-1.96*sd.est$se[i],digits=2),scientific=TRUE),";",format(signif(sd.est$EstimatedValue[i]+1.96*sd.est$se[i],digits=2),scientific=TRUE),"]")," ")}))
    rownames(sd.est) <- sd.est$parameter
    sd.est <- sd.est[rownames(to.print),-1]

    to.print <- dplyr::mutate(dplyr::mutate(sd.est,EstimatedValue = sapply(sd.est$EstimatedValue,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)})),se=sapply(sd.est$se,FUN=function(p){format(signif(p,digits=2),scientific=TRUE)}))

    print(to.print)

  }
  cat("\t",dashed.short,"\u2022 Final Scores :\n")
  cat("\t\u00B7 OFV  :",-2*remix.output$finalRes$LL,"\n")
  cat("\t\u00B7 BIC  :",BIC(remix.output),"\n")
  cat("\t\u00B7 BICc :",BICc(remix.output),"\n")
  cat("\nExecuted in",remix.output$finalRes$iter,"iterations, ",round(remix.output$finalRes$time,digits=1),"s.\n")
  cat(plain.line)
}

#' @export
summary.cvRemix <- function(cvremix.output,...){
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"

  cat("\n",plain.line)
  cat("Outputs of cv.remix algorithm from the project :\n\t",cvremix.output$info$project,"\n")
  cat(dashed.line)
  cat("Results over the grid of \u03bb:\n")
  toprint <- apply(data.frame(
                          OFV = -2*cvremix.output$LL,
                          BIC = BIC(cvremix.output),
                          BICc = BICc(cvremix.output),
                          "Positives"=sapply(cvremix.output$res,FUN=function(r){sum(r$alpha!=0)}),row.names = round(cvremix.output$lambda,digits=2)),2,FUN=function(x){round(x,digits=2)})
  print(toprint)
  cat(dashed.line)
}

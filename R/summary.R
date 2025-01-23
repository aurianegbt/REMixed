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
      if(remix.outputs$info$finalSAEM){paste0("\n\t-final SAEM computed",
      if(remix.outputs$info$test){";\n\t-final SAEM computed"})},
      ".\n")
  cat(dashed.line)
  cat("\u2022 Final Estimated Parameters :\n")
  if(is.null(remix.outputs$finalRes$standardError)){
    print(sapply(remix.outputs$finalRes$param[-which(names(remix.outputs$finalRes$param) %in% remix.outputs$info$regParam.toprint)],FUN=function(x){round(x,digits=2)}))
  }else{

  }
  cat("\t",dashed.short,"\u2022 Final Selected Biomarkers and Estimate:\n")
  if(is.null(remix.outputs$finalRes$standardError)){
    print(sapply(remix.outputs$finalRes$param[names(which(remix.outputs$finalRes$param[remix.outputs$info$regParam.toprint]!=0))],FUN=function(x){signif(x,digits=3)}))
  }else{

  }
  cat("\t",dashed.short,"\u2022 Final Scores :\n")
  cat("\t\u00B7 OFV  :",-2*remix.outputs$finalRes$LL,"\n")
  cat("\t\u00B7 BIC  :",remix.outputs$finalRes$BIC,"\n")
  cat("\t\u00B7 eBIC :",remix.outputs$finalRes$eBIC,"\n")
  cat("\nExecuted in",remix.outputs$finalRes$iter,"iterations, ",remix.outputs$finalRes$time,"s.\n")
  cat(plain.line)
}

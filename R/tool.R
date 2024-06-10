plot
# $SAEM
#   names(data)[names(data) == parameter] <- "parameter"
# changePhase <- data$iteration[which(diff(data$phase) ==
#                                       1) + 1]
# color_phase <- rgb(212, 66, 66, maxColorValue = 255)
# color_diag <- rgb(70, 130, 180, maxColorValue = 255)
# if (parameter == "convergenceIndicator") {
#   color_diag <- rgb(212, 49, 229, maxColorValue = 255)
# }
# iteration <- NULL
# p <- ggplot(data, aes(x = iteration, y = parameter)) + geom_line(color = color_diag) +
#   geom_vline(aes(xintercept = changePhase), color = color_phase)
# if (grepl("_", parameter, fixed = TRUE)) {
#   param_split <- strsplit(parameter, "_")[[1]]
#   settings$title <- parse(text = paste0(param_split[1],
#                                         "[", param_split[2], "]"))
# }
# else {
#   settings$title <- parameter
# }
# p <- lixoftTheme(p, settings)
# p <- p + xlab(NULL) + ylab(NULL)
# return(p)

# $convergence criterion

# $LL

# $

#' Get monolix demo project path.
#'
#' @return path to the monolix demo from REMix package.
#' @export
#'
#' @examples
#' print(getMLXdir())
getMLXdir <- function(){
  return(system.file("extdata", "demoMLX.mlxtran", package = "REMix"))
}

getSMLXdir <- function(){
  return(system.file("extdata", "demoSMLX.mlxtran", package = "REMix"))
}

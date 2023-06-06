#' @title Modify resources and resource availability
#' @description Implements [Resources] for the null model of resources
#' @inheritParams Resources
#' @return none
#' @export
Resources.null <- function(t, y, pars) {
  return(pars)
}

#' @title Set up parameters for the null model for resource availability
#' @param pars a [list]
#' @return none
#' @export
setup_resources_null<- function(pars){
  RApar <- list()
  class(RApar) <- 'null'
  pars$RApar <- RApar
  return(pars)
}

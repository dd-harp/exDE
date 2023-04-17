#' @title Modify parameters due to exogenous forcing
#' @description Implements [ResourceAvailability] for the null model of hydrology (do nothing)
#' @inheritParams ResourceAvailability
#' @return none
#' @export
ResourceAvailability.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model for resource availability (do nothing)
#' @param pars a [list]
#' @return none
#' @export
make_parameters_resources_null<- function(pars) {
  RApar <- list()
  class(RApar) <- 'null'
  pars$RApar <- RApar
  return(pars)
}

# specialized methods for a model with exogenous forcing

#' @title Modify parameters due to exogenous forcing
#' @description Calls functions that set up exogenous variables
#' @param t current simulation time
#' @param pars a [list]
#' @return pars a [list]
#' @export
ExogenousForcing.forced <- function(t, pars) {
  pars = Weather(t, pars)
  pars = Hydrology(t, pars)
  return(pars)
}

#' @title Make parameters for the null model of exogenous forcing (do nothing)
#' @param pars a [list]
#' @return none
#' @export
setup_exogenous_forced <- function(pars) {
  EXOpar <- list()
  class(EXOpar) <- 'forced'
  pars$EXOpar <- EXOpar
  pars = setup_weather_null(pars)
  pars = setup_hydrology_null(pars)
  return(pars)
}

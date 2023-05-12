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
  pars = ResourceAvailability(t, pars)
  return(pars)
}

#' @title Make parameters for the null model of exogenous forcing (do nothing)
#' @param pars a [list]
#' @return none
#' @export
make_parameters_exogenous_forced <- function(pars) {
  EXOpar <- list()
  class(EXOpar) <- 'forced'
  pars$EXOpar <- EXOpar
  pars$weather = make_parameters_weather_null(pars)
  pars$hydrology = make_parameters_hydrology_null(pars)
  pars$resources = make_parameters_resources_null(pars)
  return(pars)
}

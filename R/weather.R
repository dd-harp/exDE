# generic methods for exogenous forcing by weather

#' @title Methods for exogenous variables describing weather
#' @description This method dispatches on the type of `pars$WEATHER`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
Weather <- function(t, pars) {
  UseMethod("Weather", pars$WEATHER)
}

#' @title Methods for exogenous variables describing weather
#' @description Implements [Weather] for the null model (no variables)
#' @inheritParams Weather
#' @return [list]
#' @export
Weather.null <- function(t, pars) {
  return(pars)
}

#' @title Set up the null model for weather
#' @param pars a [list]
#' @return [list]
#' @export
setup_weather_null <- function(pars) {
  WEATHER <- list()
  class(WEATHER) <- 'null'
  pars$WEATHER <- WEATHER
  return(pars)
}

#' @title Methods for exogenous variables describing weather
#' @description Implements exogenous forcing by [Weather]
#' @inheritParams Weather
#' @return [list]
#' @export
Weather.forced <- function(t, pars) {
  return(pars)
}

#' @title Set up the forced model for weather
#' @param pars a [list]
#' @return [list]
#' @export
setup_weather_forced <- function(pars) {
  pars = check_abiotic(pars)
  WEATHER <- list()
  class(WEATHER) <- 'forced'
  pars$WEATHER <- WEATHER
  pars = setup_temperature_null(pars)
  pars = setup_rainfall_null(pars)
  pars = setup_relative_humidity_null(pars)
  return(pars)
}

#' @title Set up the null model for temperature
#' @param pars a [list]
#' @return [list]
#' @export
setup_temperature_null <- function(pars) {
  TEMPERATURE <- list()
  class(TEMPERATURE) <- 'null'
  pars$TEMPERATURE <- TEMPERATURE
  return(pars)
}

#' @title Set up the null model for RAINFALL
#' @param pars a [list]
#' @return [list]
#' @export
setup_rainfall_null <- function(pars) {
  RAINFALL <- list()
  class(RAINFALL) <- 'null'
  pars$RAINFALL <- RAINFALL
  return(pars)
}

#' @title Set up the null model for HUMIDITY
#' @param pars a [list]
#' @return [list]
#' @export
setup_relative_humidity_null <- function(pars) {
  HUMIDITY <- list()
  class(HUMIDITY) <- 'null'
  pars$HUMIDITY <- HUMIDITY
  return(pars)
}

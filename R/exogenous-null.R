# specialized methods for the null model of exogenous forcing

#' @title Modify parameters due to exogenous forcing
#' @description Implements [ExogenousForcing] for the null model of exogenous forcing (do nothing)
#' @inheritParams ExogenousForcing
#' @return none
#' @export
ExogenousForcing.null <- function(t, pars) {pars}

#' @title Set up the null model for exogenous forcing (do nothing)
#' @param pars a [list]
#' @return none
#' @export
setup_exogenous_null <- function(pars) {
  EXOpar <- list()
  class(EXOpar) <- 'null'
  pars$EXOpar <- EXOpar
  return(pars)
}

#' @title Modify exogenous variables describing weather
#' @description Implements [Weather] for the null model (no variables)
#' @inheritParams Weather
#' @return none
#' @export
Weather.null <- function(t, pars) {
  return(pars)
}

#' @title Set up the null model for weather (do nothing)
#' @param pars a [list]
#' @return none
#' @export
setup_weather_null <- function(pars) {
  WETpar <- list()
  class(WETpar) <- 'null'
  pars$WETpar <- WETpar
  return(pars)
}

#' @title Modify parameters due to exogenous forcing
#' @description Implements [Hydrology] for the null model of hydrology (do nothing)
#' @inheritParams Hydrology
#' @return none
#' @export
Hydrology.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model for hydrology (do nothing)
#' @param pars a [list]
#' @return none
#' @export
setup_hydrology_null <- function(pars) {
  HYpar <- list()
  class(HYpar) <- 'null'
  pars$HYpar <- HYpar
  return(pars)
}

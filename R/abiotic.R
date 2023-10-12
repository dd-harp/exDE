# generic methods for exogenous forcing by abiotic factors

#' @title Set up exogenous variables for abiotic forcing
#' @description This method dispatches on the type of `pars$ABIOTIC`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
Abiotic <- function(t, pars) {
  UseMethod("Abiotic", pars$ABIOTIC)
}

#' @title Set up exogenous variables for abiotic forcing
#' @description Implements [Abiotic] for the null model of exogenous forcing (do nothing)
#' @inheritParams Abiotic
#' @return [list]
#' @export
Abiotic.null <- function(t, pars) {pars}

#' @title Set up the null model for exogenous forcing (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_abiotic_null <- function(pars) {
  ABIOTIC <- list()
  class(ABIOTIC) <- 'null'
  pars$ABIOTIC <- ABIOTIC
  return(pars)
}

#' @title Set up exogenous variables for abiotic forcing
#' @description Implements [Abiotic] for abiotic forcing
#' @param t current simulation time
#' @param pars a [list]
#' @return pars a [list]
#' @export
Abiotic.forced <- function(t, pars) {
  pars = Weather(t, pars)
  pars = Hydrology(t, pars)
  return(pars)
}

#' @title Make parameters for the null model of abiotic forcing (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_abiotic_forced <- function(pars) {
  ABIOTIC <- list()
  class(ABIOTIC) <- 'forced'
  pars$ABIOTIC <- ABIOTIC
  pars = setup_weather_null(pars)
  pars = setup_hydrology_null(pars)
  return(pars)
}

#' @title Check abiotic
#' @param pars a [list]
#' @return [list]
#' @export
check_abiotic <- function(pars) {
  UseMethod("check_abiotic", pars$EfSz)
}

#' @title Check abiotic
#' @param pars a [list]
#' @return [list]
#' @export
check_abiotic.null <- function(pars) {
  setup_abiotic_forced(pars)
}

#' @title Check abiotic
#' @param pars a [list]
#' @return [list]
#' @export
check_abiotic.forced<- function(pars) {pars}


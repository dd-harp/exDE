# generic methods for exogenous forcing

#' @title Modify parameters due to exogenous forcing
#' @description This method dispatches on the type of `pars$EXOpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return none
#' @export
ExogenousForcing <- function(t, pars) {
  UseMethod("ExogenousForcing", pars$EXOpar)
}

#' @title Set the values of exogenous variables describing weather
#' @description This method dispatches on the type of `pars$Wpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return none
#' @export
Weather <- function(t, pars) {
  UseMethod("Weather", pars$Wpar)
}

#' @title Set the values of exogenous variables describing hydrology
#' @description This method dispatches on the type of `pars$HYpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return none
#' @export
Hydrology <- function(t, pars) {
  UseMethod("Intervene", pars$HYpar)
}




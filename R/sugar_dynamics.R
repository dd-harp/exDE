# Methods to set up variables describing exogenous forcing by sugar

#' @title Set the values of exogenous variables describing sugar
#' @description This method dispatches on the type of `pars$SUGAR`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
SugarDynamics <- function(t, pars) {
  UseMethod("SugarDynamics", pars$SUGAR)
}

#' @title Set the values of exogenous variables describing sugar
#' @description Implements [SugarDynamics] for the static model of sugar (do nothing)
#' @inheritParams SugarDynamics
#' @return [list]
#' @export
SugarDynamics.static <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the static model for sugar (do nothing)
#' @param pars a [list]
#' @param Sugar describes sugar availability
#' @return [list]
#' @export
setup_sugar_static <- function(pars, Sugar=0) {
  SUGAR <- list()
  class(SUGAR) <- 'static'
  pars$SUGAR <- SUGAR
  pars$vars$Sugar = Sugar
  return(pars)
}



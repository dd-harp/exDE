# Methods to set up variables describing exogenous forcing by sugar

#' @title Set the values of exogenous variables describing sugar
#' @description This method dispatches on the type of `pars$SUGAR`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
Sugar <- function(t, pars) {
  UseMethod("Sugar", pars$SUGAR)
}

#' @title Set the values of exogenous variables describing sugar
#' @description Implements [Sugar] for the null model of sugar (do nothing)
#' @inheritParams Sugar
#' @return [list]
#' @export
Sugar.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model for sugar (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_sugar_null <- function(pars) {
  SUGAR <- list()
  class(SUGAR) <- 'null'
  pars$SUGAR <- SUGAR
  return(pars)
}


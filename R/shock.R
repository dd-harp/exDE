# generic methods for shocks

#' @title Set up shocks
#' @description This method dispatches on the type of `pars$SHOCK`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
Shock <- function(t, pars) {
  UseMethod("Shock", pars$SHOCK)
}

#' @title Set up shocks
#' @description Implements [Shock] for the null model (do nothing)
#' @inheritParams Shock
#' @return [list]
#' @export
Shock.null <- function(t, pars) {
  return(pars)
}

#' @title Set up the null model for shocks (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_shock_null <- function(pars) {
  SHOCK <- list()
  class(SHOCK) <- 'null'
  pars$SHOCK <- SHOCK
  return(pars)
}

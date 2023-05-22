#' @title Set the values of exogenous variables describing available non-host resources
#' @description This method dispatches on the type of `pars$RApar`.
#' @param t current simulation time
#' @param y state variables
#' @param pars a [list]
#' @return none
#' @export
Resources <- function(t, y, pars) {
  UseMethod("Resources", pars$RApar)
}

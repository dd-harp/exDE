# generic methods for exogenous forcing

#' @title Modify parameters due to exogenous forcing
#' @description This method dispatches on the type of `pars$EXpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return none
#' @export
ExogenousForcing <- function(t, y, pars) {
  UseMethod("ExogenousForcing", pars$EXpar)
}

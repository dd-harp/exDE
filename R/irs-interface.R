# generic methods for IRS

#' @title Modify baseline bionomic values due to IRS
#' @description This method dispatches on the type of `pars$IRSpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
IRS <- function(t, pars) {
  UseMethod("IRS", pars$IRSpar)
}

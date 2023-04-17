#' @title Set the values of exogenous variables describing available non-host resources
#' @description This method dispatches on the type of `pars$RESpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return none
#' @export
ResourceAvailability <- function(t, pars) {
  UseMethod("ResourceAvailability", pars$RESpar)
}

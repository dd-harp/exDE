# generic methods for ATSB

#' @title Modify baseline bionomic values due to ATSB
#' @description This method dispatches on the type of `pars$ATSBpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
ATSB <- function(t, pars) {
  UseMethod("ATSB", pars$ATSBpar)
}

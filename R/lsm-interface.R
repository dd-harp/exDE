# generic methods for LSM

#' @title Modify values due to LSM
#' @description This method dispatches on the type of `pars$LSMpar`.
#' @param t current simulation time
#' @param pars an [environment]
#' @return a [list]
#' @export
LSM <- function(t, pars) {
  UseMethod("LSM", pars$LSMpar)
}

# generic methods to compute the extrinsic incubation period (EIP)

#' @title Compute the EIP
#' @description This method dispatches on the type of `pars$EIPmod`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
EIP <- function(t, pars) {
  UseMethod("EIP", pars$EIPmod)
}

#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description Implements [EIP] for the static model (the EIP is constant)
#' @inheritParams EIP
#' @return [list]
#' @export
EIP.static <- function(t, pars) {pars}

#' @title Set up the static model for control forcing (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_eip_static<- function(pars) {
  EIPmod <- list()
  class(EIPmod) <- 'static'
  pars$EIPmod <- EIPmod
  return(pars)
}

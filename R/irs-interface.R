# generic methods for IRS

#' @title Do mass house spraying (IRS)
#' @description This method dispatches on the type of `pars$IRS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
SprayHouses <- function(t, pars) {
  UseMethod("SprayHouses", pars$IRS)
}

#' @title Model the effects of IRS
#' @description This method dispatches on the type of `pars$IRS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
IRSeffects <- function(t, pars) {
  UseMethod("IRS_effects", pars$IRS)
}

#' @title Model IRS effect sizes
#' @description This method dispatches on the type of `pars$IRS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
IRSeffectSizes <- function(t, pars) {
  UseMethod("IRSeffectSizes", pars$IRS)
}

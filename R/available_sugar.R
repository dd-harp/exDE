
#' @title Update the availability of sugar
#' @description This method dispatches on the type of `pars$SUGAR`.
#' @param pars an [list]
#' @return a [list]
#' @export
AvailableSugar <- function(pars) {
  UseMethod("AvailableSugar", pars$SUGAR)
}

#' @title Compute total availability of sugar
#' @description Computes the availability of sugar for the static model (do nothing)
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableSugar.static <- function(pars){
  return(pars)
}

#' @title Compute total availability of sugar
#' @description Computes the availability of sugar
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
AvailableSugar.forced <- function(pars){
  pars$S = pars$nectar + pars$sugar_baits
  return(pars)
}

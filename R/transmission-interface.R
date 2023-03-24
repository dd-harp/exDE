# generic methods for transmission

#' @title Biting distribution matrix
#' @description This method dispatches on the type of `pars$xde`
#' @param t current simulation time
#' @param y state vector
#' @param pars, a [list]
#' @return pars, a [list]
#' @export
F_beta <- function(t, y, pars){
  UseMethod("F_beta", pars$xde)
}

#' @title Entomological inoculation rate on human strata
#' @description This method dispatches on the type of `pars$Xpar$xde`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars) {
  UseMethod("F_EIR", pars$Xpar$xde)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar$xde`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa <- function(t, y, pars) {
  UseMethod("F_kappa", pars$MYZpar$xde)
}


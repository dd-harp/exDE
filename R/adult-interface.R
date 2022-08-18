# generic methods for adult component

#' @title Entomological inoculation rate on human strata
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars) {
  UseMethod("F_EIR", pars$MYZpar)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa <- function(t, y, pars) {
  UseMethod("F_kappa", pars$MYZpar)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs <- function(t, y, pars) {
  UseMethod("F_eggs", pars$MYZpar)
}

#' @title Derivatives for adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param Lambda emergence rate of adult mosquitoes
#' @param kappa net infectiousness of human population
#' @return a [numeric] vector
#' @export
dMYZdt <- function(t, y, pars, Lambda, kappa) {
  UseMethod("dMYZdt", pars$MYZpar)
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @return the modified parameter [list]
#' @export
make_index_MYZ <- function(pars) {
  UseMethod("make_index_MYZ", pars$MYZpar)
}

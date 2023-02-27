# generic methods for demography (nested within human; \cal{H} in \cal{X})

#' @title Size of human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H <- function(t, y, pars) {
  UseMethod("F_H", pars$Hpar)
}

#' @title Size of lagged human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param lag duration of lag `t-lag`
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H_lag <- function(t, y, pars, lag) {
  UseMethod("F_H_lag", pars$Hpar)
}

#' @title Derivatives of demographic changes in human populations
#' @description This method dispatches on the type of `pars$Hpar$Hmatrix`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return see help pages for specific methods
#' @export
dHdt <- function(t, y, pars){
  UseMethod("dHdt", pars$Hpar$Hmatrix)
}

#' @title A function that computes the birth rate for human populations
#' @description This method dispatches on the type of `pars$Hpar$birthF`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return see help pages for specific methods
#' @export
Births <- function(t, y, pars){
  UseMethod("Births", pars$Hpar$birthF)
}

#' @title Add indices for human population denominators to parameter list
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @return none
#' @export
make_indices_H <- function(pars) {
  UseMethod("make_indices_H", pars$Hpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_H <- function(pars) {
  UseMethod("get_inits_H", pars$Hpar)
}

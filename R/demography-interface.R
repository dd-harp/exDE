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

#' @title Derivatives of demographic changes in human populations
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param ... additional arguments which must be of the form `D`, `Y`, etc. `D`
#' is a `nStrata` by `nStrata` demography matrix and `Y` is a state vector. There
#' should be at least one pair of these passed to the function.
#' @return see help pages for specific methods
#' @export
dHdt <- function(pars, ...) {
  UseMethod("dHdt", pars$Hpar)
}

#' @title Add indices for human population denominators to parameter list
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @return none
#' @export
make_index_H <- function(pars) {
  UseMethod("make_index_H", pars$Hpar)
}

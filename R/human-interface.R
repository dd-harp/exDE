# generic methods for human component

#' @title Entomological inoculation rate on human strata
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars) {
  UseMethod("F_EIR", pars$Xpar)
}

#' @title Size of effective infectious human population
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x <- function(t, y, pars) {
  UseMethod("F_x", pars$Xpar)
}

#' @title Size of lagged effective infectious human population
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param lag duration of lag `t-lag`
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x_lag <- function(t, y, pars, lag) {
  UseMethod("F_x_lag", pars$Xpar)
}

#' @title Derivatives for human population
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param EIR vector giving the per-capita entomological inoculation rate for each strata
#' @return a [numeric] vector
#' @export
dXdt <- function(t, y, pars, EIR) {
  UseMethod("dXdt", pars$Xpar)
}

#' @title Add indices for human population to parameter list
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return the modified parameter [list]
#' @export
make_index_X <- function(pars) {
  UseMethod("make_index_X", pars$Xpar)
}

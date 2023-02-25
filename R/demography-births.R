
#' @title The human population birth rate
#' @description This method dispatches on the type of `pars$Hpar$birthF`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return see help pages for specific methods
#' @export
F_births <- function(t, y, pars){
  UseMethod("F_births", pars$Hpar$birthF)
}

#' @title Constant human population birth rate
#' @description Implements [F_births] with a constant population birth rate
#' @inheritParams F_births
#' @return a [numeric] vector of length `nStrata`
#' @export
F_births.constant <- function(t, y, pars){
  pars$Hpar$birthrate
}

#' @title Constant per-capita human population birth rate
#' @description Implements [F_births] with a constant per-capita population birth rate
#' @inheritParams F_births
#' @return a [numeric] vector of length `nStrata`
#' @export
F_births.exp <- function(t, y, pars){
  pars$Hpar$birthrate * F_H(t,y,pars)
}

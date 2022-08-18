# generic methods for aquatic component

#' @title Number of newly emerging adults from each larval habitat
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha <- function(t, y, pars) {
  UseMethod("F_alpha", pars$Lpar)
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param eta vector giving number of eggs being laid in each larval habitat
#' @return a [numeric] vector of length `pars$L_ix`
#' @export
dLdt <- function(t, y, pars, eta) {
  UseMethod("dLdt", pars$Lpar)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$Lpar`. Adds field `L_ix`
#' to parameter list.
#' @param pars a [list]
#' @return the modified parameter [list]
#' @export
make_index_L <- function(pars) {
  UseMethod("make_index_L", pars$Lpar)
}

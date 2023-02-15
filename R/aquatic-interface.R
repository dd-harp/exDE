# generic methods for aquatic component

#' @title Number of newly emerging adults from each larval habitat
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha <- function(t, y, pars) {
  UseMethod("F_alpha", pars$Lpar)
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param eta vector giving number of eggs being laid in each larval habitat
#' @return a [numeric] vector of length `pars$L_ix`
#' @export
dLdt <- function(t, y, pars, eta) {
  UseMethod("dLdt", pars$Lpar)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$Lpar`. Adds field `L_ix`
#' to parameter list.
#' @param pars an [environment]
#' @return none
#' @export
make_indices_L <- function(pars) {
  UseMethod("make_indices_L", pars$Lpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Lpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_L <- function(pars) {
  UseMethod("get_inits_L", pars$Lpar)
}

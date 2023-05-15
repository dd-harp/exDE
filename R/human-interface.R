# generic methods for human component

#' @title Size of effective infectious human population
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X <- function(t, y, pars) {
  UseMethod("F_X", pars$Xpar)
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

#' @title A function to set up Xpar
#' @description This method dispatches on `Xname`.
#' @param pars a [list]
#' @param Xname a [character] string
#' @param Xopts a [list]
#' @return none
#' @export
setup_X = function(pars, Xname, Xopts=list()){
  class(Xname) <- Xname
  UseMethod("setup_X", Xname)
}

#' @title Add indices for human population to parameter list
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
make_indices_X <- function(pars) {
  UseMethod("make_indices_X", pars$Xpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X <- function(pars) {
  UseMethod("get_inits_X", pars$Xpar)
}

#' @title Compute the human transmitting capacity
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
HTC <- function(pars) {
  UseMethod("get_inits_X", pars$Xpar)
}


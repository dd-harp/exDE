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

#' @title Infection blocking pre-erythrocytic immunity
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b <- function(y, pars) {
  UseMethod("F_b", pars$Xpar)
}

#' @title Derivatives for human population
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param FoI vector giving the per-capita force of infection for each strata
#' @return a [numeric] vector
#' @export
dXdt <- function(t, y, pars, FoI) {
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

#' @title Parse the output of deSolve and return the variables by name in a list
#' @description This method dispatches on the type of `pars$Xpar`. Adds the variables
#' from the X model to varslist and returns it
#' @param varslist a [list] the object to be returned
#' @param deout a [matrix] of outputs from deSolve
#' @param pars a [list] that defines a model
#' @export
parse_deout_X <- function(varslist, deout, pars) {
  UseMethod("parse_deout_X", pars$Xpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X <- function(pars) {
  UseMethod("get_inits_X", pars$Xpar)
}

#' @title Set the initial values from a vector of states
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X <- function(pars, y0) {
  UseMethod("update_inits_X", pars$Xpar)
}

#' @title Compute the human transmitting capacity
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
HTC <- function(pars) {
  UseMethod("get_inits_X", pars$Xpar)
}


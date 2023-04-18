# specialized methods for a human trace model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the trace model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.trace <- function(t, y, pars) {
  with(pars$Xpar, kappa)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the trace model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.trace <- function(t, y, pars, EIR) {
  numeric(0)
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the trace model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.trace <- function(pars) {
  pars$X_ix <- integer(0)
  return(pars)
}

#' @title Make parameters for trace human model
#' @param pars an [list]
#' @param kappa net infectiousness
#' @return a [list]
#' @export
make_parameters_X_trace <- function(pars, kappa) {
  Xpar <- list()
  class(Xpar) <- c('trace')
  Xpar$kappa <- kappa
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for trace human model
#' @param pars an [environment]
#' @return none
#' @export
make_inits_X_trace <- function(pars) {
  pars$Xinits <- numeric(0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_X.trace <- function(pars){
  numeric(0)
}

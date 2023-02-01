# trace driven model of demography

#' @title Size of human population denominators
#' @description Implements [F_H] for the forced (trace) model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.trace <- function(t, y, pars) {
  pars$Hpar$H(t)
}

#' @title Size of lagged human population denominators
#' @description Implements [F_H_lag] for forced (trace) model.
#' @inheritParams F_H_lag
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H_lag.trace <- function(t, y, pars, lag) {
  pars$Hpar$H(t-lag)
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the forced (trace) model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata` or of length 0
#' @export
dHdt.trace <- function(pars, ...) {
  if (...length() > 2) {
    # being called to evaluate \dot{H}
    numeric(0)
  } else {
    # being called to evaluate a component \dot{X}
    rep(0, length(...elt(2)))
  }
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_index_H] for forced (trace) model.
#' @inheritParams make_index_H
#' @return none
#' @export
make_index_H.trace <- function(pars) {
  pars$H_ix <- integer(0)
  return(pars)
}

#' @title Make parameters for forced (trace) human demography model
#' @param pars an [environment]
#' @param H a function taking a single argument `t` and returning a vector of length
#' `nStrata`.
#' @return none
#' @export
make_parameters_demography_trace <- function(pars, H) {
  stopifnot(is.environment(pars))
  stopifnot(length(formals(H)) == 1)
  stopifnot(length(H(0)) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c('trace')
  Hpar$H <- H
  pars$Hpar <- Hpar
  return(pars)
}

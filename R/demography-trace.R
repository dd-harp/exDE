# trace driven model of demography

#' @title Size of human population denominators
#' @description Implements [F_H] for the forced (trace) model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.trace <- function(t, y, pars) {
  pars$Hpar$H(t)
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
}

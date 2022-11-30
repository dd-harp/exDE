# null model of \mathcal{H}; constant for all time

#' @title Size of human population denominators
#' @description Implements [F_H] for the null model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.null <- function(t, y, pars) {
  pars$Hpar$H
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the null model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata` or of length 0
#' @export
dHdt.null <- function(pars, ...) {
  if (...length() > 2) {
    # being called to evaluate \dot{H}
    numeric(0)
  } else {
    # being called to evaluate a component \dot{X}
    rep(0, length(...elt(2)))
  }
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_index_H] for null model.
#' @inheritParams make_index_H
#' @return none
#' @export
make_index_H.null <- function(pars) {
  pars$H_ix <- integer(0)
}

# generic methods for demography (nested within human; \cal{H} in \cal{X})

#' @title Size of human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H <- function(t, y, pars) {
  UseMethod("F_H", pars$Hpar)
}

#' @title Size of human population denominators
#' @description Implements [F_H] for the null model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.null <- function(t, y, pars) {
  pars$Hpar$H
}

#' @title Size of human population denominators
#' @description Implements [F_H] for the forced (trace) model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.trace <- function(t, y, pars) {
  pars$Hpar$H(t)
}

#' @title Size of human population denominators
#' @description Implements [F_H] for dynamic models.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.dynamic <- function(t, y, pars) {
  y[pars$H_ix]
}

#' @title Dynamics of human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param D a demography matrix, which should be of dimensions `nStrata` by `nStrata`
#' @return a [numeric] vector of length `nStrata`
#' @export
dHdt <- function(t, y, dy, pars, D = NULL) {
  UseMethod("dHdt", pars$Hpar)
}

#' @title Dynamics of human population denominators
#' @description Implements [dHdt] for the null model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata`
#' @export
dHdt.null <- function(t, y, dy, pars, D = NULL) {
  if (is.null(dy)) {
    numeric(0)
  } else {
    dy
  }
}

#' @title Dynamics of human population denominators
#' @description Implements [dHdt] for the forced (trace) model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata`
#' @export
dHdt.trace <- function(t, y, dy, pars, D = NULL) {
  if (is.null(dy)) {
    numeric(0)
  } else {
    dy
  }
}

#' @title Dynamics of human population denominators
#' @description Implements [dHdt] for dynamic models.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata`
#' @export
dHdt.dynamic <- function(t, y, dy, pars, D = NULL) {
  if (is.null(dy)) {
    D %*% y
  } else {
    dy + (D %*% y)
  }
}

#' @title Add indices for human population denominators to parameter list
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @return the modified parameter [list]
#' @export
make_index_H <- function(pars) {
  UseMethod("make_index_H", pars$Hpar)
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_index_H] for null model.
#' @inheritParams make_index_H
#' @return the modified parameter [list]
#' @export
make_index_H.null <- function(pars) {
  pars$H_ix <- integer(0)
  return(pars)
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_index_H] for forced (trace) model.
#' @inheritParams make_index_H
#' @return the modified parameter [list]
#' @export
make_index_H.trace <- function(pars) {
  pars$H_ix <- integer(0)
  return(pars)
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_index_H] for dynamic models.
#' @inheritParams make_index_H
#' @return the modified parameter [list]
#' @importFrom utils tail
#' @export
make_index_H.dynamic <- function(pars) {
  pars$H_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$H_ix, 1)
  return(pars)
}

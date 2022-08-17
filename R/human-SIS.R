# specialized methods for the human SIS model

#' @title Size of effective infectious human population
#' @description Implements [F_x] for the SIS model.
#' @inheritParams F_x
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x.SIS <- function(t, y, pars) {
  y[pars$X_ix] * pars$Xpar$c
}

#' @title Size of lagged effective infectious human population
#' @description Implements [F_x_tau] for the SIS model.
#' @inheritParams F_x_tau
#' @return a [numeric] vector of length `nStrata`
#' @importFrom deSolve lagvalue
#' @export
F_x_tau.SIS <- function(t, y, pars, tau) {
  if (t < tau) {
    X_tau <- pars$Xpar$X0
  } else {
    X_tau <- lagvalue(t = t - tau, nr = pars$X_ix)
  }
  X_tau * pars$Xpar$c
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIS <- function(t, y, pars, EIR) {
  X <- y[pars$X_ix]
  with(pars$Xpar, {
    dXdt <- diag(b*EIR) %*% (H - X) - r*X
    return(dXdt)
  })
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_index_X] for the SIS model.
#' @inheritParams make_index_X
#' @return the modified parameter [list]
#' @importFrom utils tail
#' @export
make_index_X.SIS <- function(pars) {
  pars$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$X_ix, 1)
  return(pars)
}

#' @title Make parameters for SIS human model
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list] with class `SIS`.
#' @export
make_parameters_X_SIS <- function(b, c, r) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r))
  Xpar <- list()
  class(Xpar) <- c('SIS')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  return(Xpar)
}

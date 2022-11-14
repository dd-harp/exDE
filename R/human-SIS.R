# specialized methods for the human SIS model

#' @title Entomological inoculation rate on human strata
#' @description Implements [F_EIR] for the SIS model.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR.SIS <- function(t, y, pars) {
  # Z <- y[pars$Z_ix] # may want to use wrapper compute_Z (like F_x below)
  Z <- F_Z(t, y, pars)
  f <- pars$MYZpar$f # may want to use wrapper compute_f/q
  q <- pars$MYZpar$q
  as.vector(pars$beta %*% diag(f*q, nrow = pars$nPatches) %*% Z)
}

#' @title Size of effective infectious human population
#' @description Implements [F_x] for the SIS model.
#' @inheritParams F_x
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x.SIS <- function(t, y, pars) {
  y[pars$X_ix] * pars$Xpar$c
}

#' @title Size of lagged effective infectious human population
#' @description Implements [F_x_lag] for the SIS model.
#' @inheritParams F_x_lag
#' @return a [numeric] vector of length `nStrata`
#' @importFrom deSolve lagvalue
#' @export
F_x_lag.SIS <- function(t, y, pars, lag) {
  if (t < lag) {
    X_tau <- pars$Xpar$X0
  } else {
    X_tau <- lagvalue(t = t - lag, nr = pars$X_ix)
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
    dXdt <- diag(b*EIR, nrow = pars$nStrata) %*% (H - X) - r*X
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
#' @param X0 size of infected population in each strata
#' @param H size of human population in each strata
#' @return a [list] with class `SIS`.
#' @export
make_parameters_X_SIS <- function(b, c, r, X0, H) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(X0), is.numeric(H))
  Xpar <- list()
  class(Xpar) <- c('SIS')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$X0 <- X0
  Xpar$H <- H
  return(Xpar)
}

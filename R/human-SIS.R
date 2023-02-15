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
#' @description Implements [F_x_lag] for the SIS model.
#' @inheritParams F_x_lag
#' @return a [numeric] vector of length `nStrata`
#' @importFrom deSolve lagvalue
#' @export
F_x_lag.SIS <- function(t, y, pars, lag) {
  if (t < lag) {
    X_tau <- pars$Xinits$X0
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
  H <- F_H(t, y, pars)

  with(pars$Xpar, {
    # disease dynamics
    dX <- diag(b*EIR, nrow = pars$nStrata) %*% (H - X) - r*X + dHdt(t, X, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    return(c(dX, dH))
  })
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIS model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIS <- function(pars) {
  pars$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$X_ix, 1)
  return(pars)
}

#' @title Make parameters for SIS human model
#' @param pars an [environment]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param Psi a [matrix] of dimensions `nPatches` by `nStrata`
#' @param wf vector of biting weights of length `nStrata`
#' @return none
#' @export
make_parameters_X_SIS <- function(pars, b, c, r) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r))
  Xpar <- list()
  class(Xpar) <- c('SIS')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for SIS human model
#' @param pars an [environment]
#' @param X0 size of infected population in each strata
#' @return none
#' @export
make_inits_X_SIS <- function(pars, X0) {
  stopifnot(is.numeric(X0))
  pars$Xinits <- list(X0=X0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_X.SIS <- function(pars){
  pars$Xinits$X0
}

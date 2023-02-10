# specialized methods for the human SIS model

#' @title Entomological inoculation rate on human strata
#' @description Implements [F_EIR] for the SIS model.
#' @inheritParams F_EIR
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR.SIS <- function(t, y, pars, MosyBehavior) {
  Z <- F_Z(t, y, pars)
  f <- MosyBehavior$f[1]
  q <- MosyBehavior$q[1]
  beta <- F_beta(t, y, pars)
  as.vector(beta %*% diag(f*q, nrow = pars$nPatches) %*% Z)
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

#' @title Biting distribution matrix
#' @description Implements [F_beta] for the SIS model.
#' @inheritParams F_beta
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta.SIS <- function(t, y, pars) {
  H <- F_H(t, y, pars)
  W <- as.vector(pars$Xpar$Psi %*% (pars$Xpar$wf * H))
  return(
    diag(pars$Xpar$wf, pars$nStrata) %*% t(pars$Xpar$Psi) %*% diag(1/W, pars$nPatches)
  )
}

#' @title Lagged biting distribution matrix
#' @description Implements [F_beta_lag] for the SIS model.
#' @inheritParams F_beta_lag
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta_lag.SIS <- function(t, y, pars, lag) {
  H <- F_H_lag(t, y, pars, lag)
  W <- as.vector(pars$Xpar$Psi %*% (pars$Xpar$wf * H))
  return(
    diag(pars$Xpar$wf, pars$nStrata) %*% t(pars$Xpar$Psi) %*% diag(1/W, pars$nPatches)
  )
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
#' @description Implements [make_index_X] for the SIS model.
#' @inheritParams make_index_X
#' @return none
#' @importFrom utils tail
#' @export
make_index_X.SIS <- function(pars) {
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
#' @param X0 size of infected population in each strata
#' @return none
#' @export
make_parameters_X_SIS <- function(pars, b, c, r, Psi, wf = 1, X0) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(X0))
  if (length(wf) == 1) {
    wf <- rep(wf, pars$nStrata)
  }
  stopifnot(length(wf) == pars$nStrata)
  stopifnot(nrow(Psi) == pars$nPatches)
  stopifnot(ncol(Psi) == pars$nStrata)
  Xpar <- list()
  class(Xpar) <- c('SIS')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$Psi <- Psi
  Xpar$wf <- wf
  Xpar$X0 <- X0
  pars$Xpar <- Xpar
  return(pars)
}

# specialized methods for the human SIP model

#' @title Size of effective infectious human population
#' @description Implements [F_x] for the SIP model.
#' @inheritParams F_x
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x.SIP <- function(t, y, pars) {
  y[pars$X_ix] * pars$Xpar$c
}

#' @title Size of lagged effective infectious human population
#' @description Implements [F_x_tau] for the SIP model.
#' @inheritParams F_x_tau
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x_tau.SIP <- function(t, y, pars, tau) {
  if (t < tau) {
    X_tau <- pars$Xpar$X0
  } else {
    X_tau <- deSolve::lagvalue(t = t - tau, nr = pars$X_ix)
  }
  X_tau * pars$Xpar$c
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIP model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIP <- function(t, y, pars, EIR) {
  X <- y[pars$X_ix]
  P <- y[pars$P_ix]
  with(pars$Xpar, {
    dXdt <- diag((1-rho)*b*EIR) %*% (H - X - P) - r*X
    dPdt <- diag(rho*b*EIR) %*% (H - X - P) - eta*P
    return(c(dXdt, dPdt))
  })
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_index_X] for the SIP model.
#' @inheritParams make_index_X
#' @return the modified parameter [list]
#' @importFrom utils tail
#' @export
make_index_X.SIP <- function(pars) {
  pars$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$X_ix, 1)

  pars$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$P_ix, 1)
  return(pars)
}

#' @title Make parameters for SIP human model
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param X0 size of infected population in each strata
#' @param P0 size of population protected by prophylaxis in each strata
#' @return a [list] with class `SIP`.
#' @export
make_parameters_X_SIP <- function(b, c, r, rho, eta, X0, P0) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(rho), is.numeric(eta), is.numeric(X0), is.numeric(P0))
  Xpar <- list()
  class(Xpar) <- c('SIP')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$rho <- rho
  Xpar$eta <- eta
  Xpar$X0 <- X0
  Xpar$P0 <- P0
  return(Xpar)
}
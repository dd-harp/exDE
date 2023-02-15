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
#' @description Implements [F_x_lag] for the SIP model.
#' @inheritParams F_x_lag
#' @return a [numeric] vector of length `nStrata`
#' @importFrom deSolve lagvalue
#' @export
F_x_lag.SIP <- function(t, y, pars, lag) {
  if (t < lag) {
    X_tau <- pars$Xinits$X0
  } else {
    X_tau <- lagvalue(t = t - lag, nr = pars$X_ix)
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
  H <- F_H(t, y, pars)

  with(pars$Xpar, {
    # disease dynamics
    dX <- diag((1-rho)*b*EIR, nrow = pars$nStrata) %*% (H - X - P) - r*X + dHdt(t, X, pars)
    dP <- diag(rho*b*EIR, nrow = pars$nStrata) %*% (H - X - P) - eta*P + dHdt(t, P, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    # return derivatives
    return(c(dX, dP, dH))
  })
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIP model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIP <- function(pars) {
  pars$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$X_ix, 1)

  pars$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$P_ix, 1)
  return(pars)
}

#' @title Make parameters for SIP human model
#' @param pars an [environment]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param Psi a [matrix] of dimensions `nPatches` by `nStrata`
#' @param wf vector of biting weights of length `nStrata`
#' @param X0 size of infected population in each strata
#' @param P0 size of population protected by prophylaxis in each strata
#' @return none
#' @export
make_parameters_X_SIP <- function(pars, b, c, r, rho, eta){
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(rho), is.numeric(eta))
  Xpar <- list()
  class(Xpar) <- c('SIP')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$rho <- rho
  Xpar$eta <- eta
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for SIP human model
#' @param pars an [environment]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param Psi a [matrix] of dimensions `nPatches` by `nStrata`
#' @param wf vector of biting weights of length `nStrata`
#' @param X0 size of infected population in each strata
#' @param P0 size of population protected by prophylaxis in each strata
#' @return none
#' @export
make_inits_X_SIP <- function(pars, X0, P0) {
  stopifnot(is.numeric(X0), is.numeric(P0))
  pars$Xinits = list(X0 = X0, P0 = P0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_X.SIP <- function(pars){with(pars$Xinits,{
  c(X0, P0)
})}

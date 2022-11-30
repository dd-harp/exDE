# a hybrid model tracking mean MoI for all and apparent infections

#' @title Entomological inoculation rate on human strata
#' @description Implements [F_EIR] for the hybrid MoI model.
#' @inheritParams F_EIR
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR.hMoI <- function(t, y, pars, MosyBehavior) {
  Z <- F_Z(t, y, pars)
  f <- MosyBehavior$f[1]
  q <- MosyBehavior$q[1]
  beta <- F_beta(t, y, pars)
  as.vector(beta %*% diag(f*q, nrow = pars$nPatches) %*% Z)
}

#' @title Size of effective infectious human population
#' @description Implements [F_x] for the hybrid MoI model.
#' @inheritParams F_x
#' @return a [numeric] vector of length `nStrata`
#' @importFrom stats pexp
#' @export
F_x.hMoI <- function(t, y, pars) {
  x1 <- pexp(q = y[pars$m1_ix])
  x2 <- pexp(q = y[pars$m2_ix])
  x <- (pars$Xpar$c2 * x2) + (pars$Xpar$c1 * (x1 - x2))
  return(x * as.vector(pars$Xpar$H))
}

#' @title Size of lagged effective infectious human population
#' @description Implements [F_x_lag] for the hybrid MoI model.
#' @inheritParams F_x_lag
#' @return a [numeric] vector of length `nStrata`
#' @importFrom stats pexp
#' @importFrom deSolve lagvalue
#' @export
F_x_lag.hMoI <- function(t, y, pars, lag) {
  if (t < lag) {
    m1_tau <- pars$Xpar$m10
    m2_tau <- pars$Xpar$m20
  } else {
    m1_tau <- lagvalue(t = t - lag, nr = pars$m1_ix)
    m2_tau <- lagvalue(t = t - lag, nr = pars$m2_ix)
  }
  x1_tau <- pexp(q = m1_tau)
  x2_tau <- pexp(q = m2_tau)
  x_tau <- (pars$Xpar$c2 * x2_tau) + (pars$Xpar$c1 * (x1_tau - x2_tau))
  return(x_tau * as.vector(pars$Xpar$H))
}

#' @title Biting distribution matrix
#' @description Implements [F_beta] for the hybrid MoI model.
#' @inheritParams F_beta
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta.hMoI <- function(t, y, pars) {
  H <- F_H(t, y, pars)
  W <- as.vector(pars$Xpar$Psi %*% (pars$Xpar$wf * H))
  return(
    diag(pars$Xpar$wf, pars$nStrata) %*% t(pars$Xpar$Psi) %*% diag(1/W, pars$nPatches)
  )
}

#' @title Lagged biting distribution matrix
#' @description Implements [F_beta_lag] for the hybrid MoI model.
#' @inheritParams F_beta_lag
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta_lag.hMoI <- function(t, y, pars, lag) {
  H <- F_H_lag(t, y, pars, lag)
  W <- as.vector(pars$Xpar$Psi %*% (pars$Xpar$wf * H))
  return(
    diag(pars$Xpar$wf, pars$nStrata) %*% t(pars$Xpar$Psi) %*% diag(1/W, pars$nPatches)
  )
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the hybrid MoI model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.hMoI <- function(t, y, pars, EIR) {
  m1 <- y[pars$m1_ix]
  m2 <- y[pars$m2_ix]
  with(pars$Xpar, {
    dm1dt <- b*EIR - r1*m1
    dm2dt <- b*EIR - r2*m2
    return(c(dm1dt, dm2dt))
  })
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_index_X] for the hybrid MoI model.
#' @inheritParams make_index_X
#' @return none
#' @importFrom utils tail
#' @export
make_index_X.hMoI <- function(pars) {
  pars$m1_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$m1_ix, 1)

  pars$m2_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$m2_ix, 1)
}

#' @title Make parameters for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param pars an [environment]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c1 transmission probability (efficiency) from inapparent human infections to mosquito
#' @param c2 transmission probability (efficiency) from patent human infections to mosquito
#' @param r1 recovery rate from inapparent infections
#' @param r2 recovery rate from patent infections
#' @param Psi a [matrix] of dimensions `nPatches` by `nStrata`
#' @param wf vector of biting weights of length `nStrata`
#' @param m10 mean MoI among inapparent human infections
#' @param m20 mean MoI among patent human infections
#' @return none
#' @export
make_parameters_X_hMoI <- function(pars, b, c1, c2, r1, r2, Psi, wf = 1, m10, m20) {
  stopifnot(is.numeric(b), is.numeric(c1), is.numeric(c2), is.numeric(r1), is.numeric(r2), is.numeric(m10), is.numeric(m20))
  stopifnot(is.environment(pars))
  if (length(wf) == 1) {
    wf <- rep(wf, pars$nStrata)
  }
  stopifnot(length(wf) == pars$nStrata)
  stopifnot(nrow(Psi) == pars$nPatches)
  stopifnot(ncol(Psi) == pars$nStrata)
  Xpar <- list()
  class(Xpar) <- c('hMoI')
  Xpar$b <- b
  Xpar$c1 <- c1
  Xpar$c2 <- c2
  Xpar$r1 <- r1
  Xpar$r2 <- r2
  Xpar$Psi <- Psi
  Xpar$wf <- wf
  Xpar$m10 <- m10
  Xpar$m20 <- m20
  pars$Xpar <- Xpar
}

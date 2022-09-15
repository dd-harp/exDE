# a hybrid model tracking mean MoI for all and apparent infections

#' @title Size of effective infectious human population
#' @description Implements [F_x] for the hybrid MoI model.
#' @inheritParams F_x
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x.hMoI <- function(t, y, pars) {
  x1 <- pexp(q = y[pars$m1_ix])
  x2 <- pexp(q = y[pars$m2_ix])
  x <- (pars$Xpar$c2 * x2) + (pars$Xpar$c1 * (x1 - x2))
  return(x * as.vector(pars$Xpar$H))
}

#' @title Size of lagged effective infectious human population
#' @description Implements [F_x_tau] for the hybrid MoI model.
#' @inheritParams F_x_tau
#' @return a [numeric] vector of length `nStrata`
#' @export
F_x_tau.hMoI <- function(t, y, pars, tau) {
  if (t < tau) {
    m1_tau <- pars$Xpar$m10
    m2_tau <- pars$Xpar$m20
  } else {
    m1_tau <- deSolve::lagvalue(t = t - tau, nr = pars$m1_ix)
    m2_tau <- deSolve::lagvalue(t = t - tau, nr = pars$m2_ix)
  }
  x1_tau <- pexp(q = m1_tau)
  x2_tau <- pexp(q = m2_tau)
  x_tau <- (pars$Xpar$c2 * x2_tau) + (pars$Xpar$c1 * (x1_tau - x2_tau))
  return(x_tau * as.vector(pars$Xpar$H))
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
#' @return the modified parameter [list]
#' @importFrom utils tail
#' @export
make_index_X.hMoI <- function(pars) {
  pars$m1_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$m1_ix, 1)

  pars$m2_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$m2_ix, 1)
  return(pars)
}

#' @title Make parameters for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c1 transmission probability (efficiency) from inapparent human infections to mosquito
#' @param c2 transmission probability (efficiency) from patent human infections to mosquito
#' @param r1 recovery rate from inapparent infections
#' @param r2 recovery rate from patent infections
#' @param m10 mean MoI among inapparent human infections
#' @param m20 mean MoI among patent human infections
#' @param H size of human population in each strata
#' @return a [list] with class `hMOI`.
#' @export
make_parameters_X_SIP <- function(b, c, r, rho, eta, X0, P0, H) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(rho), is.numeric(eta), is.numeric(X0), is.numeric(P0), is.numeric(H))
  Xpar <- list()
  class(Xpar) <- c('hMOI')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$rho <- rho
  Xpar$eta <- eta
  Xpar$X0 <- X0
  Xpar$P0 <- P0
  Xpar$H <- H
  return(Xpar)
}

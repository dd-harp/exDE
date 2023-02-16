# a hybrid model tracking mean MoI for all and apparent infections

#' @title Size of effective infectious human population
#' @description Implements [F_x] for the hybrid MoI model.
#' @inheritParams F_x
#' @return a [numeric] vector of length `nStrata`
#' @importFrom stats pexp
#' @export
F_x.hMoI <- function(t, y, pars) {
  H <- F_H(t, y, pars)
  x1 <- pexp(q = y[pars$m1_ix])
  x2 <- pexp(q = y[pars$m2_ix])
  x <- (pars$Xpar$c2 * x2) + (pars$Xpar$c1 * (x1 - x2))
  return(x * H)
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
    m1_tau <- pars$Xinits$m10
    m2_tau <- pars$Xinits$m20
  } else {
    m1_tau <- lagvalue(t = t - lag, nr = pars$m1_ix)
    m2_tau <- lagvalue(t = t - lag, nr = pars$m2_ix)
  }
  H <- F_H_lag(t, y, pars, lag)
  x1_tau <- pexp(q = m1_tau)
  x2_tau <- pexp(q = m2_tau)
  x_tau <- (pars$Xpar$c2 * x2_tau) + (pars$Xpar$c1 * (x1_tau - x2_tau))
  return(x_tau * H)
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
#' @description Implements [make_indices_X] for the hybrid MoI model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.hMoI <- function(pars) {
  pars$m1_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$m1_ix, 1)

  pars$m2_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$m2_ix, 1)
  return(pars)
}

#' @title Make parameters for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param pars an [environment]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c1 transmission probability (efficiency) from inapparent human infections to mosquito
#' @param c2 transmission probability (efficiency) from patent human infections to mosquito
#' @param r1 recovery rate from inapparent infections
#' @param r2 recovery rate from patent infections
#' @return none
#' @export
make_parameters_X_hMoI <- function(pars, b, c1, c2, r1, r2) {
  stopifnot(is.numeric(b), is.numeric(c1), is.numeric(c2), is.numeric(r1), is.numeric(r2))
  Xpar <- list()
  class(Xpar) <- c('hMoI')
  Xpar$b <- b
  Xpar$c1 <- c1
  Xpar$c2 <- c2
  Xpar$r1 <- r1
  Xpar$r2 <- r2
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param pars an [environment]
#' @param m10 mean MoI among inapparent human infections
#' @param m20 mean MoI among patent human infections
#' @return none
#' @export
make_inits_X_hMoI <- function(pars, m10, m20) {
  stopifnot(is.numeric(m10), is.numeric(m20))
  pars$Xinits = list(m10 = m10, m20 = m20)
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_X.hMoI <- function(pars){with(pars$Xinits,{
  c(m10, m20)
})}

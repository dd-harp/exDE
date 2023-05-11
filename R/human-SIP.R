# specialized methods for the human SIP model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIP model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIP <- function(t, y, pars) {
  with(pars$Xpar, y[X_ix] * c)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIP model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPdX <- function(t, y, pars, EIR) {

  with(pars$Xpar, {

    foi = F_foi(b*EIR, pars)
    X <- y[X_ix]
    P <- y[P_ix]
    H <- F_H(t, y, pars)

    dX <- (1-rho)*foi*(H - X - P) - (r+xi)*X
    dP <- rho*foi*(H - X - P) + xi*(H-P) - eta*P

    return(c(dX, dP))
  })
}


#' @title Compute the HTC for the SIP model
#' @description Implements [HTC] for the SIP model with demography.
#' @inheritParams HTC
#' @return a [numeric] vector
#' @export
HTC.SIP <- function(pars) {
  with(pars$Xpar,
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIP model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPdXdH <- function(t, y, pars, EIR) {

  with(pars$Xpar, {

    foi = F_foi(b*EIR, pars)
    X <- y[X_ix]
    P <- y[P_ix]
    H <- F_H(t, y, pars)

    dX <- (1-rho)*foi*(H - X - P) - (r+xi)*X + dHdt(t, X, pars)
    dP <- rho*foi*(H - X - P) + xi*(H-P) - eta*P + dHdt(t, P, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

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
  pars$Xpar$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$X_ix, 1)

  pars$Xpar$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$P_ix, 1)
  return(pars)
}

#' @title Make parameters for SIP human model
#' @param pars an [environment]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param xi background treatment rate
#' @return none
#' @export
make_parameters_X_SIP <- function(pars, b, c, r, rho, eta, xi){
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(rho), is.numeric(eta), is.numeric(xi))
  Xpar <- list()
  class(Xpar) <- c('SIP', 'SIPdX')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$rho <- rho
  Xpar$eta <- eta
  Xpar$xi <- xi
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for SIP human model
#' @param pars an [environment]
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

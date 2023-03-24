# methods for transmission involving delay differential equations

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_kappa] for when MYZ is a DDE model.
#' @inheritParams F_kappa
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa.dde <- function(t, y, pars) {
  x <- F_x(t, y, pars)
  x_lag <- F_x_lag(t, y, pars, pars$MYZpar$tau)

  kappa <- matrix(data = 0, nrow = 2, ncol = pars$nPatches)
  kappa[1, ] <- as.vector(t(pars$beta) %*% x)
  kappa[2, ] <- as.vector(t(pars$beta_lag) %*% x_lag)
  return(kappa)
}

#' @title Entomological inoculation rate on human strata
#' @description Implements [F_EIR] for when X is an DDE model.
#' @inheritParams F_EIR
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR.dde <- function(t, y, pars){
  Z <- F_Z(t, y, pars)
  Z_lag <- F_Z_lag(t, y, pars)
  f <- pars$MYZpar$f[1]
  q <- pars$MYZpar$q[1]

  EIR <- matrix(data = 0, nrow = 2, ncol = pars$nPatches)
  EIR[1, ] <-  as.vector(pars$beta %*% diag(f*q, nrow = pars$nPatches) %*% Z)
  EIR[2, ] <-  as.vector(pars$beta_lag %*% diag(f*q, nrow = pars$nPatches) %*% Z_lag)

  return(EIR)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_beta] for a DDE model.
#' @inheritParams F_beta
#' @return a [list] vector of length `nPatches`
#' @export
F_beta.dde <- function(t, y, pars) {
  pars <- make_beta(t, y, pars)
  pars <- make_beta_lag(t, y, pars, pars$MYZpar$tau)
  return(pars)
}

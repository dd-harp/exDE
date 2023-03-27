# methods for transmission involving ordinary differential equations

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_kappa] for when MYZ is an ODE model.
#' @inheritParams F_kappa
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa.ode <- function(t, y, pars) {
  x <- F_x(t, y, pars)
  kappa =  kappa_with_visitors(t(pars$beta) %*% x, pars)
  as.vector(kappa)
}

#' @title Entomological inoculation rate on human strata
#' @description Implements [F_EIR] for when X is an ODE model.
#' @inheritParams F_EIR
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR.ode <- function(t, y, pars) {
  Z <- F_Z(t, y, pars)
  f <- pars$MYZpar$f[1]
  q <- pars$MYZpar$q[1]
  as.vector(pars$beta %*% diag(f*q*loc_fqZ(pars), nrow = pars$nPatches) %*% Z)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_beta] for an ODE model.
#' @inheritParams F_beta
#' @return a [list] vector of length `nPatches`
#' @export
F_beta.ode <- function(t, y, pars) {
  pars <- make_beta(t, y, pars)
  return(pars)
}

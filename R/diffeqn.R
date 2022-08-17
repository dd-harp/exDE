# generalized spatial differential equations

#' @title Generalized spatial differential equation model
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn <- function(t, y, pars) {

  # eta: egg laying
  eggs <- F_eggs(t, y, pars)
  eta <- pars$calU %*% eggs

  # lambda: emergence of adults
  alpha <- F_alpha(t, y, pars)
  Lambda <- pars$calN %*% alpha

  # EIR: entomological inoculation rate
  EIR <- F_EIR(t, y, pars)

  # kappa: net infectiousness of humans
  kappa <- F_kappa(t, y, pars)

  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  dMYZ <- dMYZdt(t, y, pars, Lambda, kappa)
  dX <- dXdt(t, y, pars, EIR)

  return(list(c(dL, dMYZ, dX)))
}

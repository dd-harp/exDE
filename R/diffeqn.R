# generalized spatial differential equations


#' @title Generalized spatial differential equation model
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component. The arguments `EIR_delta` and `kappa_delta` are
#' for adding external forcing to the system from unmodeled sources. This can arise
#' if humans can acquire infection by traveling outside the spatial domain, and
#' arises for mosquitoes if traveling outside the spatial domain or are being infected
#' by unmodeled (non-human) sources. By default these are set to `NULL` and are
#' turned off.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param EIR_delta a vector of values to be added to the internal `EIR`
#' @param kappa_delta a vector of values to be added to the internal `kappa`
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn <- function(t, y, pars, EIR_delta = NULL, kappa_delta = NULL) {

  # weather, climate, etc
  ExogenousForcing(t, y, pars)

  # baseline mosquito feeding and mortality
  MosyBehavior0 <- MosquitoBehavior(t, y, pars)

  # mosquito feeding and mortality under control
  MosyBehavior <- VectorControl(t, y, pars, MosyBehavior0)

  # eta: egg laying
  eggs <- F_eggs(t, y, pars)
  eta <- pars$calU %*% eggs

  # lambda: emergence of adults
  alpha <- F_alpha(t, y, pars)
  Lambda <- pars$calN %*% alpha

  # bloodmeal
  pars <- F_beta(t, y, pars)

  # EIR: entomological inoculation rate
  EIR <- F_EIR(t, y, pars, MosyBehavior)
  if (!is.null(EIR_delta)) {
    EIR <- EIR + EIR_delta
  }

  # kappa: net infectiousness of humans
  kappa <- F_kappa(t, y, pars)
  if (!is.null(kappa_delta)) {
    kappa <- kappa + kappa_delta
  }

  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  dMYZ <- dMYZdt(t, y, pars, Lambda, kappa, MosyBehavior)
  dX <- dXdt(t, y, pars, EIR)

  return(list(c(dL, dMYZ, dX)))
}

#' @title Generalized spatial differential equation model (mosquito only)
#' @description Mirrors [exDE::xDE_diffeqn] but only includes the adult and aquatic
#' mosquito components.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param kappa a vector
#' @param MosyBehavior a [list] emulating the output of [exDE::MosquitoBehavior] for
#' the appropriate adult mosquito model
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_mosy <- function(t, y, pars, kappa, MosyBehavior) {

  # weather, climate, etc
  ExogenousForcing(t, y, pars)

  # baseline mosquito feeding and mortality
  MosyBehavior0 <- MosquitoBehavior(t, y, pars)

  # mosquito feeding and mortality under control
  MosyBehavior <- VectorControl(t, y, pars, MosyBehavior0)

  # eta: egg laying
  eggs <- F_eggs(t, y, pars)
  eta <- pars$calU %*% eggs

  # lambda: emergence of adults
  alpha <- F_alpha(t, y, pars)
  Lambda <- pars$calN %*% alpha

  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  dMYZ <- dMYZdt(t, y, pars, Lambda, kappa, MosyBehavior)

  return(list(c(dL, dMYZ)))
}

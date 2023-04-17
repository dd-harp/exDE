# generalized spatial differential equations


#' @title Generalized spatial differential equation model
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn <- function(t, y, pars) {

  # set the values of exogenous forcing variables
  pars <- ExogenousForcing(t, pars)

  # set baseline mosquito bionomic parameters
  pars <- MosquitoBehavior(t, y, pars)

  # modify baseline mosquito bionomic parameters
  pars <- VectorControl(t, y, pars)

  # importation from outside the spatial domain
  pars <- Import(t, y, pars)

  # eta: egg laying
  eggs <- F_eggs(t, y, pars)
  eta <- pars$calU %*% eggs

  # lambda: emergence of adults
  alpha <- F_alpha(t, y, pars)
  Lambda <- pars$calN %*% alpha

  # blood feeding & mixing
  beta <- F_beta(t, y, pars)

  # EIR: entomological inoculation rate
  EIR <- F_EIR(t, y, pars, beta)

  # kappa: net infectiousness of humans
  kappa <- F_kappa(t, y, pars, beta)


  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  dMYZ <- dMYZdt(t, y, pars, Lambda, kappa)
  dX <- dXdt(t, y, pars, EIR)

  return(list(c(dL, dMYZ, dX)))
}

#' @title Generalized spatial differential equation model (mosquito only)
#' @description Mirrors [exDE::xDE_diffeqn] but only includes the adult and aquatic
#' mosquito components.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param kappa a vector or NULL object
#' the appropriate adult mosquito model
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_mosy <- function(t, y, pars, kappa=NULL) {

  # set the values of exogenous forcing variables
  pars <- ExogenousForcing(t, pars)

  # set baseline mosquito bionomic parameters
  pars <- MosquitoBehavior(t, y, pars)

  # modify baseline mosquito bionomic parameters
  pars <- VectorControl(t, y, pars)

  # eta: egg laying
  eggs <- F_eggs(t, y, pars)
  eta <- pars$calU %*% eggs

  # lambda: emergence of adults
  alpha <- F_alpha(t, y, pars)
  Lambda <- pars$calN %*% alpha

  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  dM <- dMYZdt(t, y, pars, Lambda, kappa)
  return(list(c(dL, dM)))
}

#' @title Differential equation models for human cohorts
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param a age of a cohort
#' @param y state vector
#' @param pars an [environment]
#' @param F_eir a trace function that returns the eir
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_cohort <- function(a, y, pars, F_eir) {

  # EIR: entomological inoculation rate trace
  eir <- F_eir(a, pars)

  # state derivatives
  dX <- dXdt(a, y, pars, eir)

  return(list(c(dX)))
}

#' @title Differential equation models for aquatic
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param F_eta trace function that returns eggs laid
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_aquatic <- function(t, y, pars, F_eta) {

  # set the values of exogenous forcing variables
  pars <- ExogenousForcing(t, pars)

  # modify baseline mosquito bionomic parameters
  pars <- LSM(t, pars)

  eta <- F_eta(t, pars)

  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  return(list(c(dL)))
}

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

  # set the values of exogenous forcing variables
  pars <- Abiotic(t, pars)
  pars <- Shock(t, pars)
  pars <- Control(t, y, pars)
  pars <- Behavior(t, y, pars)
  pars <- Visitors(t, pars)
  pars <- VectorControlEffects(t, y, pars)
  pars <- Resources(t, y, pars)

  # set and modify the baseline bionomic parameters
  pars <- MBionomics(t, y, pars, 1)
  pars <- LBionomics(t, y, pars, 1)
  pars <- VectorControlEffectSizes(t, y, pars)

  # eta: egg laying
  eta <- LayEggs(t, y, pars, 1)

  # Lambda: emergence of adults
  Lambda <- pars$calN %*% F_alpha(t, y, pars, 1)

  # blood feeding & mixing
  beta <- F_beta(t, y, pars, 1)

  # EIR: entomological inoculation rate
  EIR <- F_EIR(t, y, pars, beta, 1)

  # FoI: force of infection
  FoI <- Exposure(t, y, pars, EIR, 1)

  # kappa: net infectiousness of humans
  kappa <- F_kappa(t, y, pars, beta, 1)

  # state derivatives
  dL <- dLdt(t, y, pars, eta, 1)
  dMYZ <- dMYZdt(t, y, pars, Lambda, kappa, 1)
  dX <- dXdt(t, y, pars, FoI, 1)

  return(list(c(dL, dMYZ, dX)))
}

#' @title Differential equations isolating the humans, forced with Ztrace
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_human <- function(t, y, pars) {

  # set the values of exogenous forcing variables
  pars <- Abiotic(t, pars)
  pars <- Shock(t,  pars)
  pars <- Control(t, y, pars)
  pars <- Behavior(t, y, pars)
  pars <- Resources(t, y, pars)

  # set and modify the baseline mosquito bionomic parameters
  pars <- MBionomics(t, y, pars, 1)
  pars <- VectorControlEffectSizes(t, y, pars)

  # blood feeding & mixing
  beta <- F_beta(t, y, pars, 1)

  # EIR: entomological inoculation rate
  EIR <- F_EIR(t, y, pars, beta, 1)

  # FoI: force of infection
  FoI <- Exposure(t, y, pars, EIR, 1)

  # state derivatives
  dX <- dXdt(t, y, pars, FoI, 1)

  return(list(c(dX)))
}



#' @title Generalized spatial differential equation model (mosquito only)
#' @description Mirrors [exDE::xDE_diffeqn] but only includes the adult and aquatic
#' mosquito components.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' the appropriate adult mosquito model
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_mosy <- function(t, y, pars) {

  # set the values of exogenous forcing variables
  pars <- Abiotic(t, pars)
  pars <- Shock(t, pars)
  pars <- Control(t, y, pars)
  pars <- Behavior(t, y, pars)
  pars <- Resources(t, y, pars)

  # set baseline mosquito bionomic parameters
  pars <- MBionomics(t, y, pars, 1)
  pars <- LBionomics(t, y, pars, 1)
  pars <- VectorControlEffectSizes(t, y, pars)

  # eta: egg laying
  eta <- LayEggs(t, y, pars, 1)

  # Lambda: emergence of adults
  Lambda <- pars$calN %*% F_alpha(t, y, pars, 1)

  kappa = pars$kappa
  # state derivatives
  dL <- dLdt(t, y, pars, eta, 1)
  dM <- dMYZdt(t, y, pars, Lambda, kappa, 1)

  return(list(c(dL, dM)))
}

#' @title Differential equation models for human cohorts
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param a age of a cohort
#' @param y state vector
#' @param pars a [list]
#' @param F_eir a trace function that returns the eir
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_cohort <- function(a, y, pars, F_eir) {

  # EIR: entomological inoculation rate trace
  EIR <- F_eir(a, pars)*pars$Hpar[[1]]$wts_f

  # FoI: force of infection
  FoI <- Exposure(a, y, pars, EIR, 1)

  # state derivatives
  dX <- dXdt(a, y, pars, FoI, 1)

  return(list(c(dX)))
}

#' @title Differential equation models for aquatic mosquito populations
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_aquatic <- function(t, y, pars) {

  # set the values of exogenous forcing variables
  pars <- Abiotic(t, pars)
  pars <- Shock(t, pars)
  pars <- Control(t, y, pars)
  pars <- Resources(t, y, pars)

  # modify baseline mosquito bionomic parameters
  pars <- LBionomics(t, y, pars,1)
  pars <- VectorControlEffectSizes(t, y, pars)

  # egg laying
  eta <- LayEggs(t, y, pars, 1)

  # state derivatives
  dL <- dLdt(t, y, pars, eta, 1)
  return(list(c(dL)))
}

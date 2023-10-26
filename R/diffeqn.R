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
  pars <- EIP(t, pars)
  pars <- MBionomics(t, y, pars)
  pars <- LBionomics(t, y, pars)
  pars <- VectorControlEffectSizes(t, y, pars)

  # eta: egg laying
  eggs <- F_eggs(t, y, pars)
  eta  <- pars$calU %*% eggs

  # lambda: emergence of adults
  alpha  <- F_alpha(t, y, pars)
  Lambda <- pars$calN %*% alpha

  # blood feeding & mixing
  beta <- F_beta(t, y, pars)

  # EIR: entomological inoculation rate
  EIR <- F_EIR(t, y, pars, beta)

  # FoI: force of infection
  FoI <- Exposure(t, y, pars, EIR)

  # kappa: net infectiousness of humans
  kappa <- F_kappa(t, y, pars, beta)

  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  dMYZ <- dMYZdt(t, y, pars, Lambda, kappa)
  dX <- dXdt(t, y, pars, FoI)

  if(pars$eir_out==TRUE) {EIRt = EIR}else{EIRt = numeric(0)}
  if(pars$fqZ_out==TRUE) {fqZt = F_fqZ(t, y, pars)}else{fqZt = numeric(0)}
  if(pars$NI_out==TRUE) {NIt = F_X(t, y, pars)/F_H(t, y, pars)}else{NIt = numeric(0)}
  if(pars$kappa_out==TRUE) {kappat = kappa}else{kappat = numeric(0)}
  return(list(c(dL, dMYZ, dX, EIRt, fqZt, NIt, kappat)))
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
  pars <- EIP(t, pars)

  # set and modify the baseline mosquito bionomic parameters
  pars <- EIP(t, pars)
  pars <- MBionomics(t, y, pars)
  pars <- VectorControlEffectSizes(t, y, pars)

  # blood feeding & mixing
  beta <- F_beta(t, y, pars)

  # EIR: entomological inoculation rate
  EIR <- F_EIR(t, y, pars, beta)

  # FoI: force of infection
  FoI <- Exposure(t, y, pars, EIR)

  # state derivatives
  dX <- dXdt(t, y, pars, FoI)

  if(pars$eir_out==TRUE) {EIRt = EIR}else{EIRt = numeric(0)}
  if(pars$fqZ_out==TRUE) {fqZt = F_fqZ(t, y, pars)}else{fqZt = numeric(0)}
  if(pars$NI_out==TRUE) {NIt = F_X(t, y, pars)/F_H(t, y, pars)}else{NIt = numeric(0)}
  if(pars$kappa_out==TRUE) {kappat = kappa}else{kappat = numeric(0)}
  return(list(c(dX, EIRt, fqZt, NIt, kappat)))
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
  pars <- MBionomics(t, y, pars)
  pars <- LBionomics(t, y, pars)
  pars <- VectorControlEffectSizes(t, y, pars)

  # eta: egg laying
  eggs <- F_eggs(t, y, pars)
  eta <- pars$calU %*% eggs

  # lambda: emergence of adults
  alpha <- F_alpha(t, y, pars)
  Lambda <- pars$calN %*% alpha

  kappa = pars$kappa
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
#' @param pars a [list]
#' @param F_eir a trace function that returns the eir
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_cohort <- function(a, y, pars, F_eir) {

  # EIR: entomological inoculation rate trace
  EIR <- F_eir(a, pars)

  # FoI: force of infection
  FoI <- Exposure(a, y, pars, EIR)

  # state derivatives
  dX <- dXdt(a, y, pars, FoI)

  if(pars$eir_out==TRUE) {EIRt = EIR}else{EIRt = numeric(0)}
  if(pars$NI_out==TRUE) {NIt = F_X(a, y, pars)/F_H(a, y, pars)}else{NIt = numeric(0)}
  return(list(c(dX, EIRt, NIt)))
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
  pars <- LBionomics(t, y, pars)
  pars <- VectorControlEffectSizes(t, y, pars)

  # egg laying
  eta <- F_eggs(t, y, pars)

  # state derivatives
  dL <- dLdt(t, y, pars, eta)
  return(list(c(dL)))
}

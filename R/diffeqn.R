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
  pars <- Bionomics(t, y, pars)
  pars <- VectorControlEffectSizes(t, y, pars)

  # egg laying: compute eta
  pars <- EggLaying(t, y, pars)

  # emergence: compute Lambda
  pars <- Emergence(t, y, pars)

  # compute beta, EIR, and kappa
  pars <- Transmission(t, y, pars)

  # compute the FoI
  pars <- Exposure(t, y, pars)

  # compute derivatives
  dL <- dLdt(t, y, pars, 1)
  dMYZ <- dMYZdt(t, y, pars, 1)
  if(pars$nVectors > 1)
    for(s in 2:pars$nVectors){
      dL <- c(dL, dLdt(t, y, pars, s))
      dMYZ <- c(dMYZ, dMYZdt(t, y, pars, s))
    }

  dX <- dXdt(t, y, pars, 1)
  if(pars$nHosts > 1)
    for(i in 2:pars$nHosts)
      dX <- c(dX, dXdt(t, y, pars, i))


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

  # compute beta, EIR, and kappa
  pars <- Transmission(t, y, pars)

  # compute the FoI
  pars <- Exposure(t, y, pars)

  # state derivatives
  dX <- dXdt(t, y, pars, 1)
  if(pars$nHosts > 1)
    for(i in 2:pars$nHosts)
      dX <- c(dX, dXdt(t, y, pars, i))

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
  #pars <- Resources(t, y, pars)

  # set baseline mosquito bionomic parameters
  pars <- Bionomics(t, y, pars)
  pars <- VectorControlEffectSizes(t, y, pars)

  # egg laying: compute eta
  pars <- EggLaying(t, y, pars)

  # emergence: compute Lambda
  pars <- Emergence(t, y, pars)

  # state derivatives
  dL <- dLdt(t, y, pars, 1)
  dM <- dMYZdt(t, y, pars, 1)
  if (pars$nVectors > 1)
    for(s in 2:pars$nVectors){
      dL <- c(dL, dLdt(t, y, pars, s))
      dM <- c(dM, dMYZdt(t, y, pars, s))
    }

  return(list(c(dL, dM)))
}

#' @title Differential equation models for human cohorts
#' @description Compute derivatives for [deSolve::ode] or [deSolve::dede] using
#' generic methods for each model component.
#' @param a age of a cohort
#' @param y state vector
#' @param pars a [list]
#' @param F_eir a trace function that returns the eir as a function of time
#' @return a [list] containing the vector of all state derivatives
#' @export
xDE_diffeqn_cohort <- function(a, y, pars, F_eir) {

  # EIR: entomological inoculation rate trace
  pars$EIR[[1]] <- with(pars$EIRpar, F_eir(a, bday, scale))*pars$BFpar$relativeBitingRate[[1]][[1]]

  # FoI: force of infection
  pars <- Exposure(a, y, pars)

  # state derivatives
  dX <- dXdt(t, y, pars, 1)
  if(pars$nHosts > 1)
    for(i in 2:pars$nHosts)
      dX <- c(dX, dXdt(t, y, pars, i))

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
  pars <- HabitatDynamics(t, pars)

  # modify baseline mosquito bionomic parameters
  pars <- LBionomics(t, y, pars, 1)
  pars <- VectorControlEffectSizes(t, y, pars)

  # egg laying: compute eta

  pars$eggs_laid[[1]] = F_eggs(t, y, pars, 1)

  # state derivatives
  dL <- dLdt(t, y, pars, 1)
  if(pars$nVectors > 1)
    for(s in 1:pars$nVectors)
      dL <- c(dL, dLdt(t, y, pars, s))

  return(list(c(dL)))
}

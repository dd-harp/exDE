
#' @title Compute other variables at time t
#' @description Compute everything but the derivatives for the generalized
#' spatial differential equation model
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return **pars** a [list]
#' @export
compute_vars_full <- function(t, y, pars) {

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

  return(pars)
}

#' @title Compute other variables at time t
#' @description Compute everything but the derivatives for a human-only
#' differential equation model
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return **pars** a [list]
#' @export
compute_vars_human <- function(t, y, pars) {

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

  return(pars)
}

#' @title Compute other variables at time t
#' @description Compute everything but the derivatives for a mosquito-only
#' differential equation model
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' the appropriate adult mosquito model
#' @return **pars** a [list]
#' @export
compute_vars_mosy <- function(t, y, pars) {

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

  return(pars)
}

#' @title Differential equation models for human cohorts
#' @description Compute everything but the derivatives for a human cohort
#' differential equation model
#' @param a age of a cohort
#' @param y state vector
#' @param pars a [list]
#' @param F_eir a trace function that returns the eir
#' @return **pars** a [list]
#' @export
compute_vars_cohort <- function(a, y, pars, F_eir) {

  # EIR: entomological inoculation rate trace
  pars$EIR[[1]] <- F_eir(a, pars)*pars$BFpar$relativeBitingRate[[1]][[1]]

  # FoI: force of infection
  pars <- Exposure(a, y, pars)

  return(pars)
}

#' @title Differential equation models for aquatic mosquito populations
#' @description Compute everything but the derivatives for an aquatic mosquito-only
#' differential equation model
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return **pars** a [list]
#' @export
compute_vars_aqua <- function(t, y, pars) {

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

  return(pars)
}

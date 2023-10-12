# specialized methods for models with dynamically changing resource availability

#' @title Methods for resources
#' @description Implements [Resources]
#' @inheritParams Resources
#' @return [list]
#' @export
Resources.forced <- function(t, y, pars) {

  # Availability of Blood Hosts
  pars = OtherBloodHosts(t, pars)
  pars = HumanAvailability(t, y, pars)

  # Availability of Blood Hosts
  pars = HabitatDynamics(t, pars)
  pars = HabitatAvailability(pars)

  # Availability of Blood Hosts
  pars = Sugar(t, pars)
  pars = SugarAvailability(pars)

  pars$local_frac = compute_local_frac(pars)

  return(pars)
}


#' @title Compute availability of local humans for blood feeding
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
HumanAvailability.forced <- function(t, y, pars){
  H = F_H(t, y, pars)
  pars$W = with(pars$Hpar, TaR %*% (wts_f*H))
  return(pars)
}

#' @title Compute total availability of aquatic habitats
#' @description Computes the availability of aquatic habitats
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
HabitatAvailability.forced <- function(pars){
  habitats = with(pars, calN %*% searchQ)
  pars$Q =  habitats + pars$ovitraps + pars$non_habitats
  pars$eggs_laid = habitats/pars$Q
  pars$calU = make_calU(pars$calN, pars$searchQ)
  return(pars)
}

#' @title Compute total availability of sugar
#' @description Computes the availability of sugar
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
SugarAvailability.forced <- function(pars){
  pars$S = pars$nectar + pars$sugar_baits
  return(pars)
}

#' @title Set up a model for mass medical
#' @param pars a [list]
#' @return [list]
#' @export
setup_resources_forced <- function(pars) {
  pars = setup_control(pars)
  RESOURCES <- list()
  class(RESOURCES) <- 'forced'
  pars$RESOURCES <- RESOURCES
  pars = setup_habitat_dynamics_null(pars)
  pars = setup_sugar_null(pars)
  pars = setup_other_blood_hosts_null(pars)
  return(pars)
}

#' @title Compute the local fraction
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param pars a [list]
#' @return pars a [list]
#' @export
compute_local_frac <- function(pars){
  pars$local_frac = with(pars, W/(W+Visitors))
  return(pars)
}

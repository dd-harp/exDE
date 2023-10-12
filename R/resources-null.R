#' @title Modify resources and resource availability
#' @description Implements [Resources] for the null model of resources
#' @inheritParams Resources
#' @return none
#' @export
Resources.null <- function(t, y, pars) {
  return(pars)
}

#' @title Compute availability of local humans for blood feeding
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
HumanAvailability.null <- function(t, y, pars){
  return(pars)
}

#' @title Compute total availability of aquatic habitats
#' @description Computes the availability of aquatic habitats for the null model (do nothing)
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
HabitatAvailability.null <- function(pars){
  return(pars)
}

#' @title Compute total availability of sugar
#' @description Computes the availability of sugar for the null model (do nothing)
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
SugarAvailability.null <- function(pars){
  return(pars)
}

#' @title Set up parameters for the null model for resource availability
#' @param pars a [list]
#' @return none
#' @export
setup_resources_null<- function(pars){
  RESOURCES <- list()
  class(RESOURCES) <- 'null'
  pars$RESOURCES <- RESOURCES
  return(pars)
}

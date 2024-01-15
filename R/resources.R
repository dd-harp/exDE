# Methods to compute resource dynamics and availability

#' @title Set the values of exogenous variables describing available mosquito resources
#' @description This method dispatches on the type of `pars$RESOURCES`.
#' @param t current simulation time
#' @param y vector of state variables
#' @param pars a [list]
#' @return none
#' @export
Resources <- function(t, y, pars) {
  UseMethod("Resources", pars$RESOURCES)
}

#' @title Modify resources and resource availability
#' @description Implements [Resources] for the null model of resources
#' @inheritParams Resources
#' @return none
#' @export
Resources.static <- function(t, y, pars) {
  return(pars)
}

#' @title Modify resources and resource availability
#' @description Implements [Resources] for the null model of resources
#' @inheritParams Resources
#' @return none
#' @export
Resources.setup <- function(t, y, pars) {

  pars = OtherBloodHosts(t, pars)
  pars = AvailableBlood(t, y, pars)

  pars = HabitatDynamics(t, pars)
  pars = AvailableHabitat(pars)

  pars = SugarDynamics(t, pars)
  pars = AvailableSugar(pars)

  class(pars$RESOURCES) <- "static"

  return(pars)
}

#' @title Methods for resources
#' @description Implements [Resources]
#' @inheritParams Resources
#' @return [list]
#' @export
Resources.forced <- function(t, y, pars) {

  pars = OtherBloodHosts(t, pars)
  pars = AvailableBlood(t, y, pars)

  pars = HabitatDynamics(t, pars)
  pars = AvailableHabitat(pars)

  pars = SugarDynamics(t, pars)
  pars = AvailableSugar(pars)

  return(pars)
}

#' @title Set up parameters for the null model for resource availability
#' @param pars a [list]
#' @return none
#' @export
setup_resources_null <- function(pars){
  RESOURCES <- list()
  class(RESOURCES) <- 'static'
  pars$RESOURCES <- RESOURCES
  pars <- setup_other_blood_hosts_static(pars)
  return(pars)
}

#' @title Set up parameters for the null model for resource availability
#' @param pars a [list]
#' @return none
#' @export
setup_resources_static <- function(pars){
  RESOURCES <- list()
  class(RESOURCES) <- 'setup'
  pars$RESOURCES <- RESOURCES

  pars$vars$Q = list()
  pars$vars$non_habitats = 0
  pars$vars$ovitraps = 0

  pars = setup_sugar_static(pars)
  pars = setup_other_blood_hosts_static(pars)
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
  pars = setup_sugar_static(pars)
  pars = setup_other_blood_hosts_static(pars)
  return(pars)
}

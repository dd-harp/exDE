# Methods to set up variables describing exogenous forcing by habitat_dynamics

#' @title Habitat Dynamics and Searching
#' @description Set the values of habitat search weights and other exogenous variables describing habitat_dynamics. This method dispatches on the type of `pars$HABITAT_DYNAMICS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
HabitatDynamics <- function(t, pars) {
  UseMethod("HabitatDynamics", pars$HABITAT_DYNAMICS)
}

#' @title Set the values of habitat search weights and other exogenous variables describing habitat_dynamics
#' @description Implements [HabitatDynamics] for the static model of habitat_dynamics (do nothing)
#' @inheritParams HabitatDynamics
#' @return [list]
#' @export
HabitatDynamics.static <- function(t, pars) {
  return(pars)
}

#' @title Setup the egg laying object
#' @description Sets up the egg-deposition matrix calU for the s^th species
#' @param pars the model object
#' @return a [list] vector
#' @export
setup_habitat_dynamics_static = function(pars){
  up <- list()
  class(up) <- "static"
  pars$HABITAT_DYNAMICS <- up
  return(pars)
}

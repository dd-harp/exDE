# Methods to set up variables describing exogenous forcing by habitat_dynamics

#' @title Set the values of exogenous variables describing habitat_dynamics
#' @description This method dispatches on the type of `pars$HABITAT_DYNAMICS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
HabitatDynamics <- function(t, pars) {
  UseMethod("HabitatDynamics", pars$HABITAT_DYNAMICS)
}

#' @title Set the values of exogenous variables describing habitat dynamics
#' @description Implements [HabitatDynamics] for the null model of habitat_dynamics (do nothing)
#' @inheritParams HabitatDynamics
#' @return [list]
#' @export
HabitatDynamics.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model for habitat dynamics (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_habitat_dynamics_null <- function(pars) {
  HABITAT_DYNAMICS <- list()
  class(HABITAT_DYNAMICS) <- 'null'
  pars$HABITAT_DYNAMICS <- HABITAT_DYNAMICS
  return(pars)
}

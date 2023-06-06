#' @title Set the values of exogenous variables describing available mosquito resources
#' @description This method dispatches on the type of `pars$RApar`.
#' @param t current simulation time
#' @param y state variables
#' @param pars a [list]
#' @return none
#' @export
Resources <- function(t, y, pars) {
  UseMethod("Resources", pars$RApar)
}

#' @title Update the availability of blood hosts
#' @description This method dispatches on the type of `pars$B`. It
#' computes availability of all blood hosts at time `t`
#' @param t current simulation time
#' @param y state variables
#' @param pars an [list]
#' @return a [list]
#' @export
update_BloodFeeding <- function(t, y, pars) {
  UseMethod("update_BloodFeeding", pars$B)
}

#' @title Update the availability of habitats and egg laying
#' @description This method dispatches on the type of `pars$Q`. It
#' computes availability of all blood hosts at time `t`
#' @param t current simulation time
#' @param y state variables
#' @param pars an [list]
#' @return a [list]
#' @export
update_EggLaying <- function(t, y, pars) {
  UseMethod("update_EggLaying", pars$Q)
}

#' @title Update host availability
#' @description This method dispatches on the type of `pars$W`. It should
#' update the availability of the pathogen's hosts for blood feeding at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
update_W <- function(t, y, pars) {
  UseMethod("update_W", pars$W)
}

#' @title Update the availability of other blood hosts
#' @description This method dispatches on the type of `pars$O`. It
#' computes availability of blood hosts at time `t`. Other blood hosts
#' include any vertebrate host that does not serve as a host for the pathogen.
#' @param t current simulation time
#' @param pars an [list]
#' @return a [list]
#' @export
update_O <- function(t, pars) {
  UseMethod("update_O", pars$O)
}

#' @title Update the total availability of sugar
#' @description This method dispatches on the type of `pars$S`. It should
#' update, and if necessary compute, the availability of sugar at time `t`
#' @param t current simulation time
#' @param pars an [list]
#' @return a [list]
#' @export
update_S <- function(t, pars) {
  UseMethod("update_S", pars$S)
}


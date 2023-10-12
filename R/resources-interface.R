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

#' @title Update the availability of blood hosts
#' @description This method dispatches on the type of `pars$RESOURCES`.
#' @param t current simulation time
#' @param y vector of state variables
#' @param pars an [list]
#' @return a [list]
#' @export
HumanAvailability <- function(t, y, pars) {
  UseMethod("HumanAvailability", pars$RESOURCES)
}

#' @title Update the availability of aquatic habitats
#' @description This method dispatches on the type of `pars$RESOURCES`.
#' @param pars an [list]
#' @return a [list]
#' @export
HabitatAvailability <- function(pars) {
  UseMethod("HabitatAvailability", pars$RESOURCES)
}

#' @title Update the availability of sugar
#' @description This method dispatches on the type of `pars$RESOURCES`.
#' @param pars an [list]
#' @return a [list]
#' @export
SugarAvailability <- function(pars) {
  UseMethod("SugarAvailability", pars$RESOURCES)
}


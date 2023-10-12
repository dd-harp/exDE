# generic methods for distributing interventions at clinics

#' @title Methods for distributing interventions during clinical visits
#' @description This method dispatches on the type of `pars$CLINIC`.
#' @param t current simulation time
#' @param y state variables
#' @param pars a [list]
#' @return [list]
#' @export
Clinic <- function(t, y, pars) {
  UseMethod("Clinic", pars$CLINIC)
}

#' @title Methods for distributing interventions during clinical visits
#' @description Implements [Clinic] for the null model (do nothing)
#' @inheritParams Clinic
#' @return [list]
#' @export
Clinic.null <- function(t, y, pars) {
  return(pars)
}

#' @title Set up the null model for clinical distribution (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_clinic_null <- function(pars) {
  CLINIC <- list()
  class(CLINIC) <- 'null'
  pars$CLINIC <- CLINIC
  return(pars)
}


#' @title Set up a model for clinic
#' @param pars a [list]
#' @return [list]
#' @export
setup_clinic_forced <- function(pars) {
  pars = setup_control(pars)
  CLINIC <- list()
  class(CLINIC) <- 'forced'
  pars$CLINIC <- CLINIC
  return(pars)
}

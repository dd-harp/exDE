# generic methods for all kinds of disease control measures

#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description This method dispatches on the type of `pars$CONTROL`.
#' @param t current simulation time
#' @param y state variables
#' @param pars a [list]
#' @return [list]
#' @export
Control <- function(t, y, pars) {
  UseMethod("Control", pars$CONTROL)
}

#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description Implements [Control] for the null model (do nothing)
#' @inheritParams Control
#' @return [list]
#' @export
Control.null <- function(t, y, pars) {pars}

#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description Implements [Control] for the static model; after setting up, do nothing
#' @inheritParams Control
#' @return [list]
#' @export
Control.static <- function(t, y, pars) {pars}

#' @title Set up the null model for control forcing (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_control_null <- function(pars) {
  CONTROL <- list()
  class(CONTROL) <- 'null'
  pars$CONTROL <- CONTROL
  return(pars)
}

#' @title Setup CONTROL with forcing
#' @param pars a [list]
#' @return [list]
#' @export
setup_control <- function(pars){
  UseMethod("setup_control", pars$CONTROL)
}

#' @title Setup control
#' @param pars a [list]
#' @return [list]
#' @export
setup_control.null <- function(pars) {
  setup_control_forced(pars)
}

#' @title Setup control
#' @param pars a [list]
#' @return [list]
#' @export
setup_control.forced<- function(pars) {pars}

#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description Implements [Control] for a model with some control
#' @param t current simulation time
#' @param y vector of state variables
#' @param pars a [list]
#' @return pars a [list]
#' @export
Control.forced <- function(t, y, pars) {
  pars = MassMedical(t, y, pars)
  pars = Clinic(t, y, pars)
  pars = VectorControl(t, y, pars)
  return(pars)
}

#' @title Set up a model with some control
#' @param pars a [list]
#' @return [list]
#' @export
setup_control_forced <- function(pars) {
  CONTROL <- list()
  class(CONTROL) <- 'forced'
  pars$CONTROL <- CONTROL
  pars = setup_mass_medical_null(pars)
  pars = setup_clinic_null(pars)
  pars = setup_vc_null(pars)
  return(pars)
}


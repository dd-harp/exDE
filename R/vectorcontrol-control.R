# specialized methods for the control model of vector control

#' @title Modify baseline values due to vector control
#' @description Implements [VectorControl] for the control model of vector control (do nothing)
#' @inheritParams VectorControl
#' @return a named [list]
#' @export
VectorControl.control <- function(t, y, pars) {
  pars = BedNets(t, pars)
  pars = IRS(t, pars)
  pars = ATSB(t, pars)
  pars = LSM(t, pars)
  return(pars)
}

#' @title Make parameters for the control model of vector control (do nothing)
#' @param pars an [environment]
#' @return none
#' @export
make_parameters_vc_control <- function(pars) {
  VCpar <- list()
  class(VCpar) <- 'control'
  pars$VCpar <- VCpar
  pars <- make_parameters_itn_null(pars)
  pars <- make_parameters_irs_null(pars)
  pars <- make_parameters_atsb_null(pars)
  pars <- make_parameters_lsm_null(pars)
  return(pars)
}

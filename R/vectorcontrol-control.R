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
#' @param pars a [list]
#' @return none
#' @export
setup_vc_control <- function(pars) {
  VCpar <- list()
  class(VCpar) <- 'control'
  pars$VCpar <- VCpar
  pars <- setup_itn_null(pars)
  pars <- setup_irs_null(pars)
  pars <- setup_atsb_null(pars)
  pars <- setup_lsm_null(pars)
  return(pars)
}

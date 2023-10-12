# specialized methods for the null model of LSM

#' @title Set up the null model for LSM
#' @description This method dispatches on the type of `pars$LSM`
#' @inheritParams TreatHabitats
#' @return a [list]
#' @export
TreatHabitats.null <- function(t, pars) {
  pars
}

#' @title Modify effects of LSM, the null model
#' @description This method dispatches on the type of `pars$LSM`
#' @inheritParams LSMeffects
#' @return a [list]
#' @export
LSMeffects.null <- function(t, y, pars) {
   pars
}

#' @title Modify effects of LSM, the null model
#' @description This method dispatches on the type of `pars$LSM`
#' @inheritParams LSMeffectSizes
#' @return a [list]
#' @export
LSMeffectSizes.null <- function(t, y, pars) {
   pars
}

#' @title Make parameters for the null model of LSM (do nothing)
#' @param pars a [list]
#' @return a [list]
#' @export
setup_lsm_null <- function(pars) {
  LSM <- list()
  class(LSM) <- 'null'
  return(pars)
}

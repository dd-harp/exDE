# specialized methods for the null model of LSM

#' @title Modify baseline values due to LSM
#' @description Implements [LSM] for the null model of LSM (do nothing)
#' @inheritParams LSM
#' @return a [list]
#' @export
LSM.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model of LSM (do nothing)
#' @param pars a [list]
#' @return a [list]
#' @export
make_parameters_lsm_null <- function(pars) {
  LSMpar <- list()
  class(LSMpar) <- 'null'
  pars$LSMpar <- LSMpar
  return(pars)
}

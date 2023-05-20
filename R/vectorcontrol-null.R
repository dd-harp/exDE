# specialized methods for the null model of vector control

#' @title Modify baseline values due to vector control
#' @description Implements [VectorControl] for the null model of vector control (do nothing)
#' @inheritParams VectorControl
#' @return a named [list]
#' @export
VectorControl.null <- function(t, y, pars) {
  return(pars)
}

#' @title Set up the vector control null model (do nothing)
#' @param pars a [list]
#' @return none
#' @export
setup_vc_null <- function(pars) {
  VCpar <- list()
  class(VCpar) <- 'null'
  pars$VCpar <- VCpar
  return(pars)
}

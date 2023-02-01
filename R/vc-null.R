# specialized methods for the null model of vector control

#' @title Modify baseline values due to vector control
#' @description Implements [VectorControl] for the null model of vector control (do nothing)
#' @inheritParams VectorControl
#' @return a named [list]
#' @export
VectorControl.null <- function(t, y, pars, MosyBehavior) {
  return(MosyBehavior)
}

#' @title Make parameters for the null model of vector control (do nothing)
#' @param pars an [environment]
#' @return none
#' @export
make_parameters_vc_null <- function(pars) {
  stopifnot(is.environment(pars))
  VCpar <- list()
  class(VCpar) <- 'null'
  pars$VCpar <- VCpar
  return(pars)
}

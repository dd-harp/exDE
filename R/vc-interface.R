# generic methods for vector control component

#' @title Modify baseline values due to vector control
#' @description This method dispatches on the type of `pars$VCpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [list]
#' @export
VectorControl <- function(t, y, pars) {
  UseMethod("VectorControl", pars$VCpar)
}

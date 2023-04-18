# generic methods for vector control component

#' @title Modify baseline values due to vector control
#' @description This method dispatches on the type of `pars$VCpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @param y state vector
#' @return a [list]
#' @export
VectorControl <- function(t, y, pars) {
  UseMethod("VectorControl", pars$VCpar)
}


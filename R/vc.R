# generic methods for vector control component

#' @title Modify baseline values due to vector control
#' @description This method dispatches on the type of `pars$VCpar`. It takes the
#' baseline `MosyBehavior` values and modifies them, potentially at multiple time
#' points for models with delay.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param MosyBehavior values returned by [exDE::MosquitoBehavior]
#' @return a [list]
#' @export
VectorControl <- function(t, y, pars, MosyBehavior) {
  UseMethod("VectorControl", pars$VCpar)
}

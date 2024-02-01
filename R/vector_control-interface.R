# generic methods for vector control component

#' @title Distribute vector control
#' @description This method dispatches on the type of `pars$VECTOR_CONTROL`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
VectorControl <- function(t, y, pars) {
  UseMethod("VectorControl", pars$VECTOR_CONTROL)
}

#' @title Vector control: durability & effects
#' @description This method dispatches on the type of `pars$VECTOR_CONTROL`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
VectorControlEffects <- function(t, y, pars) {
  UseMethod("VectorControlEffects", pars$VECTOR_CONTROL)
}

#' @title Vector control effect sizes
#' @description This method dispatches on the type of `pars$VECTOR_CONTROL`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
VectorControlEffectSizes <- function(t, y, pars) {
  UseMethod("VectorControlEffectSizes", pars$VECTOR_CONTROL)
}

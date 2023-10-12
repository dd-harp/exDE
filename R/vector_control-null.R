# specialized methods for the null model of vector control

#' @title Distribute vector control
#' @description Implements [VectorControl] for the null model of vector control (do nothing)
#' @inheritParams VectorControl
#' @return a named [list]
#' @export
VectorControl.null <- function(t, y, pars) {
  return(pars)
}

#' @title Vector control: durability & effects
#' @description Implements [VectorControlEffects] for the null model of vector control (do nothing)
#' @inheritParams VectorControlEffects
#' @return a named [list]
#' @export
VectorControlEffects.null <- function(t, y, pars) {
  return(pars)
}

#' @title Vector control effect sizes
#' @description Implements [VectorControlEffectSizes] for the null model of vector control (do nothing)
#' @inheritParams VectorControlEffectSizes
#' @return a named [list]
#' @export
VectorControlEffectSizes.null <- function(t, y, pars) {
  return(pars)
}

#' @title Distribute vector control, the null model
#' @param pars a [list]
#' @return none
#' @export
setup_vc_null <- function(pars) {
  VECTOR_CONTROL <- list()
  class(VECTOR_CONTROL) <- 'null'
  pars$VECTOR_CONTROL <- VECTOR_CONTROL
  return(pars)
}

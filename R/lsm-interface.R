# generic methods for LSM

#' @title Modify values due to treat habitats as part of LSM,
#' called by Control->VectorControl
#' @description This method dispatches on the type of `pars$LSM`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
TreatHabitats <- function(t, pars) {
  UseMethod("TreatHabitats", pars$LSM)
}

#' @title Modify effects of LSM
#' @description This method dispatches on the type of `pars$LSM`.
#' @param t current simulation time
#' @param y the state of the system
#' @param pars a [list]
#' @return a [list]
#' @export
LSM_Effects <- function(t, pars) {
  UseMethod("LSM_Effects", pars$LSM)
}

#' @title Compute effect sizes of LSM
#' @description This method dispatches on the type of `pars$LSM`.
#' @param t current simulation time
#' @param y the state of the system
#' @param pars a [list]
#' @return a [list]
#' @export
LSM_EffectSizes <- function(t, pars) {
  UseMethod("LSM_EffectSizes", pars$LSM)
}

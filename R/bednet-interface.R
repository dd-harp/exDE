# generic methods for bed nets

#' @title Distribute bed nets, called from Control(VectorControl)
#' @description This method dispatches on the type of `pars$ITNdist`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
DistributeBedNets<- function(t, pars) {
  UseMethod("DistributeBedNets", pars$ITNdist)
}

#' @title Model bed net loss, called from Control(VectorControl)
#' @description This method dispatches on the type of `pars$ITNown`.
#' @param t current simulation time
#' @param y vector of state variables
#' @param pars a [list]
#' @return a [list]
#' @export
OwnBedNet <- function(t, y, pars) {
  UseMethod("OwnBedNet", pars$ITNown)
}

#' @title Model bed net usage, called from Behavior
#' @description This method dispatches on the type of `pars$ITNuse`.
#' @param t current simulation time
#' @param y vector of state variables
#' @param pars a [list]
#' @return a [list]
#' @export
UseBedNet <- function(t, y, pars) {
  UseMethod("UseBedNets", pars$ITNuse)
}

#' @title Modify variables or parameters, called from VectorControlEffects
#' @description This method dispatches on the type of `pars$ITNeff`.
#' @param t current simulation time
#' @param pars a [list]
#' @param s the vector species index
#' @return a [list]
#' @export
BedNetEffects <- function(t, pars, s) {
  UseMethod("BedNetEffects", pars$ITNeff)
}

#' @title Modify baseline bionomic parameters, called from VectorControlEffectSizes
#' @description This method dispatches on the type of `pars$ITNefsz`.
#' @param t current simulation time
#' @param pars a [list]
#' @param s the vector species index
#' @return a [list]
#' @export
BedNetEffectSizes <- function(t, pars, s) {
  UseMethod("BedNetEffectSizes", pars$ITNefsz)
}

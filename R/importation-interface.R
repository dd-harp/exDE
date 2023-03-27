# generic methods for malaria importation

#' @title Malaria importation through infected visitors, etc.
#' @description This method dispatches on the type of `pars$MIpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return pars a [list]
#' @export
MalariaImportation <- function(t, pars) {
  UseMethod("MalariaImportation", pars$MIpar)
}

#' @title Compute kappa with visitors
#' @description This method dispatches on the type of `pars$MIpar`.
#' @param kappa, computed without visitors
#' @param pars a [list]
#' @return none
#' @export
kappa_with_visitors <- function(kappa, pars) {
  UseMethod("kappa_with_visitors", pars$MIpar)
}

#' @title Compute the local fraction, with visitors
#' @description This method dispatches on the type of `pars$MIpar`.
#' @param pars a [list]
#' @return none
#' @export
loc_fqZ <- function(pars) {
  UseMethod("loc_fqZ", pars$MIpar)
}

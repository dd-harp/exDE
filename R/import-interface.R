# generic methods for parasite / pathogen importation

#' @title Importation through visitors
#' @description This method dispatches on the type of `pars$Ipar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return pars a [list]
#' @export
Import <- function(t, y, pars) {
  UseMethod("Import", pars$Ipar)
}

#' @title Correct kappa for models with visitors
#' @description This method dispatches on the type of `pars$Ipar`.
#' @param kappa the net infectiousness among local strata
#' @param pars a [list]
#' @return pars a [list]
#' @export
kappa_local <- function(kappa, pars) {
  UseMethod("kappa_local", pars$Ipar)
}

#' @title Correct fqZ for models with visitors
#' @description This method dispatches on the type of `pars$Ipar`.
#' @param fqZ blood meals per patch, per day
#' @param pars a [list]
#' @return pars a [list]
#' @export
fqZ_local <- function(fqZ, pars) {
  UseMethod("fqZ_local", pars$Ipar)
}


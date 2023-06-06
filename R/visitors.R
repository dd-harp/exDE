# generic methods for parasite / pathogen importation by visitors

#' @title Visitors
#' @description This method dispatches on the type of `pars$Ipar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return pars a [list]
#' @export
Visitors <- function(t, pars) {
  UseMethod("Visitors", pars$Ipar)
}


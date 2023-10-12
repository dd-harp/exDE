# generic methods for oviposition traps

#' @title Methods for oviposition traps
#' @description This method dispatches on the type of `pars$OVITRAPS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
OviTraps <- function(t, pars) {
  UseMethod("OviTraps", pars$OVITRAPS)
}

#' @title Methods for oviposition traps
#' @description Implements [OviTraps] for the null model (do nothing)
#' @inheritParams OviTraps
#' @return [list]
#' @export
OviTraps.null <- function(t, pars) {
  return(pars)
}

#' @title Set up the null model for oviposition traps (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_ovitraps_null <- function(pars) {
  OVITRAPS <- list()
  class(OVITRAPS) <- 'null'
  pars$OVITRAPS <- OVITRAPS
  return(pars)
}

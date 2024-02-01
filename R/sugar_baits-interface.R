# methods for sugar baits

#' @title Methods for distributing sugar baits
#' @description This method dispatches on the type of `pars$SUGAR_BAITS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
SugarBaits <- function(t, pars) {
  UseMethod("SugarBaits", pars$SUGAR_BAITS)
}

#' @title Methods for the durability and effects of the sugar baits
#' @description This method dispatches on the type of `pars$SUGAR_BAITS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
SugarBaitEffects <- function(t, pars) {
  UseMethod("SugarBaitEffects", pars$SUGAR_BAITS)
}

#' @title Methods for the effect sizes of the sugar baits
#' @description This method dispatches on the type of `pars$SUGAR_BAITS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
SugarBaitEffectSizes <- function(t, pars) {
  UseMethod("SugarBaitEffectSizes", pars$SUGAR_BAITS)
}

#' @title Set up the null model for sugar baits (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_sugar_baits_null <- function(pars) {
  SUGAR_BAITS <- list()
  class(SUGAR_BAITS) <- 'null'
  pars$SUGAR_BAITS <- SUGAR_BAITS
  return(pars)
}

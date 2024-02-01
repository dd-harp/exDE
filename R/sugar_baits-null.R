
#' @title Methods for sugar baits
#' @description Implements [SugarBaits] for the null model (do nothing)
#' @inheritParams SugarBaits
#' @return [list]
#' @export
SugarBaits.null <- function(t, pars) {
  return(pars)
}

#' @title Methods for the effects of the sugar baits
#' @description Implements [SugarBaitEffects] for the null model (do nothing)
#' @inheritParams SugarBaitEffects
#' @return [list]
#' @export
SugarBaitEffects.null <- function(t, pars) {
  return(pars)
}

#' @title Methods for the effect sizes of the sugar baits
#' @description This method dispatches on the type of `pars$SUGAR_BAITS`.
#' @inheritParams SugarBaitEffectSizes
#' @return [list]
#' @export
SugarBaitEffectSizes.null <- function(t, pars) {
  return(pars)
}

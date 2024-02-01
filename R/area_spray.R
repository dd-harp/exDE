# generic methods for mass spraying

#' @title Methods for mass spraying
#' @description This method dispatches on the type of `pars$AREA_SPRAY`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
AreaSpray <- function(t, pars) {
  UseMethod("AreaSpray", pars$AREA_SPRAY)
}

#' @title Methods for mass spraying
#' @description This method dispatches on the type of `pars$AREA_SPRAY`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
AreaSprayEffects <- function(t, pars) {
  UseMethod("AreaSprayEffects", pars$AREA_SPRAY)
}

#' @title Methods for mass spraying
#' @description This method dispatches on the type of `pars$AREA_SPRAY`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
AreaSprayEffectSizes <- function(t, pars) {
  UseMethod("AreaSprayEffectSizes", pars$AREA_SPRAY)
}


#' @title Set up mass spraying
#' @description Implements [AreaSpray] for the null model (do nothing)
#' @inheritParams AreaSpray
#' @return [list]
#' @export
AreaSpray.null <- function(t, pars) {
  return(pars)
}

#' @title Set up mass spraying
#' @description Implements [AreaSprayEffects] for the null model (do nothing)
#' @inheritParams AreaSprayEffects
#' @return [list]
#' @export
AreaSprayEffects.null <- function(t, pars) {
  return(pars)
}

#' @title Set up mass spraying
#' @description Implements [AreaSprayEffectSizes] for the null model (do nothing)
#' @inheritParams AreaSprayEffectSizes
#' @return [list]
#' @export
AreaSprayEffectSizes.null <- function(t, pars) {
  return(pars)
}

#' @title Set up the null model for area spraying (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_area_spray_null <- function(pars) {
  AREA_SPRAY <- list()
  class(AREA_SPRAY) <- 'null'
  pars$AREA_SPRAY <- AREA_SPRAY
  return(pars)
}

# specialized methods for the null model of bed nets

#' @title Distribute bed nets
#' @description Implements [DistributeBedNets] for the null model of bed nets (do nothing)
#' @inheritParams DistributeBedNets
#' @return a [list]
#' @export
DistributeBedNets.null <- function(t, pars) {
  pars
}

#' @title Bed net ownership
#' @description Implements [OwnBedNet] for the null model of bed nets (do nothing)
#' @inheritParams OwnBedNet
#' @return a [list]
#' @export
OwnBedNet.null <- function(t, y, pars) {
  pars
}

#' @title Bed net ownership
#' @description Implements [UseBedNet] for the null model of bed nets (do nothing)
#' @inheritParams UseBedNet
#' @return a [list]
#' @export
UseBedNet.null <- function(t, y, pars) {
  pars
}

#' @title Bed net ownership
#' @description Implements [BedNetEffects] for the null model of bed nets (do nothing)
#' @inheritParams BedNetEffects
#' @return a [list]
#' @export
BedNetEffects.null <- function(t, pars, s) {
  pars
}

#' @title Bed net ownership
#' @description Implements [BedNetEffectSizes] for the null model of bed nets (do nothing)
#' @inheritParams BedNetEffectSizes
#' @return a [list]
#' @export
BedNetEffectSizes.null <- function(t, pars,s) {
  pars
}

#' @title Make parameters for the null model of bed nets (do nothing)
#' @param pars a [list]
#' @return a [list]
#' @export
setup_itn_null <- function(pars) {
  ITN<- list()
  class(ITN)   <- 'null'
  pars$ITNdist <- ITN
  pars$ITNown  <- ITN
  pars$ITNuse  <- ITN
  pars$ITNeff  <- ITN
  pars$ITNefsz <- ITN
  return(pars)
}

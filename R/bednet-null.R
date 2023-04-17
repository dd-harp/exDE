# specialized methods for the null model of bed nets

#' @title Modify baseline values due to bed nets
#' @description Implements [BedNets] for the null model of bed nets (do nothing)
#' @inheritParams BedNets
#' @return a [list]
#' @export
BedNets.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model of bed nets (do nothing)
#' @param pars a [list]
#' @return a [list]
#' @export
make_parameters_itn_null <- function(pars) {
  ITNpar <- list()
  class(ITNpar) <- 'null'
  pars$ITNpar <- ITNpar
  return(pars)
}

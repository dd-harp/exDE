# generic methods for vector control component

#' @title Modify baseline bionomic values due to bed nets
#' @description This method dispatches on the type of `pars$ITNpar`.
#' @param t current simulation time
#' @param pars an [environment]
#' @return a [list]
#' @export
BedNets <- function(t, pars) {
  UseMethod("BedNets", pars$ITNpar)
}

# generic methods for a model of travel malaria

#' @title Simulate travel malaria
#' @description This method dispatches on the type of `pars$TRAVEL`.
#' @param t current simulation time
#' @param pars a [list]
#' @return the travel FoI, a [numeric] value
#' @export
travel_malaria <- function(t, pars) {
  UseMethod("travel_malaria", pars$TRAVEL)
}


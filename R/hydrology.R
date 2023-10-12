# Methods to set up variables describing exogenous forcing by hydrology

#' @title Set the values of exogenous variables describing hydrology
#' @description This method dispatches on the type of `pars$HYDROLOGY`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
Hydrology <- function(t, pars) {
  UseMethod("Hydrology", pars$HYDROLOGY)
}

#' @title Set the values of exogenous variables describing hydrology
#' @description Implements [Hydrology] for the null model of hydrology (do nothing)
#' @inheritParams Hydrology
#' @return [list]
#' @export
Hydrology.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model for hydrology (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_hydrology_null <- function(pars) {
  HYDROLOGY <- list()
  class(HYDROLOGY) <- 'null'
  pars$HYDROLOGY <- HYDROLOGY
  return(pars)
}


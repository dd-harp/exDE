# specialized methods for the null model of IRS
# generic methods for IRS

#' @title Do mass house spraying (IRS)
#' @description This method dispatches on the type of `pars$IRS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
SprayHouses.null <- function(t, pars) {pars}

#' @title Model the effects of IRS
#' @description This method dispatches on the type of `pars$IRS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
IRSeffects.null <- function(t, pars){pars}

#' @title Model IRS effect sizes
#' @description This method dispatches on the type of `pars$IRS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
IRSeffectSizes.null <- function(t, pars){pars}

#' @title Make parameters for the null model of IRS (do nothing)
#' @param pars a [list]
#' @return a [list]
#' @export
setup_irs_null <- function(pars) {
  IRS <- list()
  class(IRS) <- 'null'
  pars$IRS <- IRS
  return(pars)
}

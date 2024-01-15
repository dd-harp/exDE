# Methods to set up variables describing exogenous forcing by other blood hosts

#' @title Set the values of exogenous variables describing other blood hosts
#' @description This method dispatches on the type of `pars$OTHER_BLOOD`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
OtherBloodHosts <- function(t, pars) {
  UseMethod("OtherBloodHosts", pars$OTHER_BLOOD)
}

#' @title Set the values of exogenous variables describing other blood hosts
#' @description Implements [OtherBloodHosts] for the static model of other_blood_hosts (do nothing)
#' @inheritParams OtherBloodHosts
#' @return [list]
#' @export
OtherBloodHosts.static <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the static model for other blood hosts (do nothing)
#' @param pars a [list]
#' @param Other the availability of other blood hosts
#' @return [list]
#' @export
setup_other_blood_hosts_static <- function(pars, Other=0) {
  OTHER_BLOOD <- list()
  class(OTHER_BLOOD) <- 'static'
  pars$vars$Other[[1]] = Other
  pars$OTHER_BLOOD <- OTHER_BLOOD
  return(pars)
}


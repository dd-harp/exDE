# specialized methods for the null model of IRS

#' @title Modify baseline values due to IRS
#' @description Implements [IRS] for the null model of IRS (do nothing)
#' @inheritParams IRS
#' @return a [list]
#' @export
IRS.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model of IRS (do nothing)
#' @param pars a [list]
#' @return a [list]
#' @export
make_parameters_irs_null <- function(pars) {
  IRSpar <- list()
  class(IRSpar) <- 'null'
  pars$IRSpar <- IRSpar
  return(pars)
}

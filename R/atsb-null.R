# specialized methods for the null model of ATSB

#' @title Modify baseline values due to ATSB
#' @description Implements [ATSB] for the null model of ATSB (do nothing)
#' @inheritParams ATSB
#' @return a [list]
#' @export
ATSB.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model of ATSB (do nothing)
#' @param pars a [list]
#' @return a [list]
#' @export
setup_atsb_null <- function(pars) {
  ATSBpar <- list()
  class(ATSBpar) <- 'null'
  pars$ATSBpar <- ATSBpar
  return(pars)
}

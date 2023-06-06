# specialized methods for the visitors null

#' @title Visitors, a null model
#' @description Implements [Visitors] for the null model (do nothing)
#' @inheritParams Visitors
#' @return a named [list]
#' @export
Visitors.null <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the null model visitors (no visitors)
#' @param pars a [list]
#' @return none
#' @export
setup_visitors_null <- function(pars) {

  Ipar <- list()
  class(Ipar) <- 'null'
  pars$Ipar <- Ipar

  pars$local_frac = 1
  pars$Visitors = 0
  pars$x_visitors = 0

  return(pars)
}

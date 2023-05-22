# specialized methods for the null model of parasite / pathogen importation

#' @title Visitors, a null model
#' @description Implements [Visitors] for the null model (do nothing)
#' @inheritParams Visitors
#' @return a named [list]
#' @export
Visitors.null <- function(t, pars) {
  return(pars)
}

#' @title Parasite / pathogen importation, the null model
#' @description Implements [kappa_local] for the null model (do nothing)
#' @inheritParams kappa_local
#' @return kappa, a numeric [vector]
#' @export
kappa_local.null<- function(kappa, pars) {
  return(kappa)
}

#' @title Parasite / pathogen importation, the null model
#' @description Implements [fqZ_local] for the null model (do nothing)
#' @inheritParams fqZ_local
#' @return fqZ, a numeric [vector]
#' @export
fqZ_local.null<- function(fqZ, pars) {
  return(fqZ)
}

#' @title Make parameters for the null model visitors (no visitors)
#' @param pars a [list]
#' @return none
#' @export
setup_visitors_null <- function(pars) {
  Ipar <- list()
  class(Ipar) <- 'null'
  pars$Ipar <- Ipar
  pars$Visitors=0
  return(pars)
}

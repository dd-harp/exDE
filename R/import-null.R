# specialized methods for the null model of parasite / pathogen importation

#' @title Visitors, a null model
#' @description Implements [Visitors] for the null model (do nothing)
#' @inheritParams Visitors
#' @return a named [list]
#' @export
Visitors.null <- function(t, pars) {
  return(pars)
}

#' @title travel_foi, a null model
#' @description Implements [travel_foi] for the null model (do nothing)
#' @inheritParams travel_foi
#' @return a [numeric]
#' @export
travel_foi.null <- function(t, pars) {
  return(0)
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

#' @title A function to set up malaria importation
#' @description Implements [setup_MalariaImportation] for the null model
#' @inheritParams setup_MalariaImportation
#' @return a [list]
#' @export
setup_MalariaImportation.null = function(pars, IMname, IMopts =list()){
  Ipar <- list()
  class(Ipar) <- 'null'
  pars$Ipar <- Ipar
  pars$Visitors=0
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
  pars$Visitors=0
  return(pars)
}

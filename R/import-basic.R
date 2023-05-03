# specialized methods for the basic model of parasite / pathogen importation

#' @title Importation, the basic model
#' @description Implements [Import] for the basic model of importation (do nothing)
#' @inheritParams Import
#' @return a named [list]
#' @export
Import.basic <- function(t, y, pars) {
  return(pars)
}

#' @title Parasite / pathogen importation, the basic model
#' @description Implements [kappa_local] for the basic model (do nothing)
#' @inheritParams kappa_local
#' @return kappa, a numeric [vector]
#' @export
kappa_local.basic<- function(kappa, pars) {
  with(pars, local_frac*kappa + (1-local_frac)*x_visitors)
  return(kappa)
}

#' @title Parasite / pathogen importation, the basic model
#' @description Implements [fqZ_local] for the basic model (do nothing)
#' @inheritParams fqZ_local
#' @return fqZ, a numeric [vector]
#' @export
fqZ_local.basic<- function(fqZ, pars) {
  return(pars$local_frac*fqZ)
}

#' @title Make parameters for the basic model for parasite / pathogen importation
#' @param pars a [list]
#' @param local_frac the fraction of bites taken on residents
#' @param x_visitors the probability a bite on a visitor infects a mosquito
#' @param Visitors is the availability of visitors for blood feeding
#' @return none
#' @export
make_parameters_import_basic <- function(pars, local_frac=1, x_visitors=0, Visitors=0) {
  Ipar <- list()
  class(Ipar) <- 'basic'
  pars$Ipar <- Ipar
  pars$x_visitors = x_visitors
  pars$local_frac = local_frac
  pars$Visitors = 0
  return(pars)
}

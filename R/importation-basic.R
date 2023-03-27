# specialized methods for the basic model of malaria importation

#' @title Malaria importation, the basic model
#' @description Implements [MalariaImportation] for the basic model of vector control (do nothing)
#' @inheritParams MalariaImportation
#' @return a named [list]
#' @export
MalariaImportation.basic <- function(t, pars) {
  return(pars)
}

#' @title kappa with visitors, the basic model
#' @description Implements [kappa_with_visitors] for the basic model of vector control (do nothing)
#' @inheritParams kappa_with_visitors
#' @return a numeric [vector]
#' @export
kappa_with_visitors.basic <- function(kappa, pars) {
  with(pars$MIpars,{
     return(local_fraction*kappa + (1-local_fraction)*x_delta)
})}

#' @title fqZ with visitors, the basic model
#' @description Implements [loc_fqZ] for the basic model of vector control (do nothing)
#' @inheritParams loc_fqZ
#' @return a numeric [vector]
#' @export
loc_fqZ.basic <- function(pars) {
  return(pars$MIpar$local_fraction)
}

#' @title Make parameters for the basic model for malaria importation
#' @param pars a [list]
#' @param local_fraction a numeric [vector] describing the fraction local
#' @param x_delta a numeric [vector] describing the density of infectious visitors
#' @return none
#' @export
make_parameters_mi_basic <- function(pars, local_fraction, x_delta) {
  MIpar <- list()
  class(MIpar) <- 'basic'
  MIpar$local_fraction = local_fraction
  MIpar$x_delta = x_delta
  pars$MIpar <- MIpar
  return(pars)
}

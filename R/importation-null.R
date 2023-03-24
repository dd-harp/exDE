# specialized methods for the null model of malaria importation

#' @title Malaria importation, the null model
#' @description Implements [MalariaImportation] for the null model of vector control (do nothing)
#' @inheritParams MalariaImportation
#' @return a named [list]
#' @export
MalariaImportation.null <- function(t, pars) {
  return(pars)
}

#' @title kappa with visitors, the null model
#' @description Implements [kappa_with_visitors] for the null model of vector control (do nothing)
#' @inheritParams kappa_with_visitors
#' @return a numeric [vector]
#' @export
kappa_with_visitors.null <- function(kappa, pars) {
  return(kappa)
}

#' @title fqZ with visitors, the null model
#' @description Implements [loc_fqZ] for the null model of vector control (do nothing)
#' @inheritParams loc_fqZ
#' @return a numeric [vector]
#' @export
loc_fqZ.null <- function(pars) {1}

#' @title Make parameters for the null model for malaria importation
#' @param pars a [list]
#' @return none
#' @export
make_parameters_mi_null <- function(pars) {
  MIpar <- list()
  class(MIpar) <- 'null'
  pars$MIpar <- MIpar
  return(pars)
}

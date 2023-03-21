# methods for transmission involving ordinary differential equations

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_kappa] for when MYZ is an trace model.
#' @inheritParams F_kappa
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa.trace <- function(t, y, pars) {
  numeric(0)
}

#' @title Entomological inoculation rate on human strata
#' @description Implements [F_EIR] for when X is an trace model.
#' @inheritParams F_EIR
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR.trace <- function(t, y, pars, MosyBehavior){
  numeric(0)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_beta] for an trace model.
#' @inheritParams F_beta
#' @return a [list] vector of length `nPatches`
#' @export
F_beta.trace <- function(t, y, pars) {
  numeric(0)
}

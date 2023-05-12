# generic methods for demography (nested within human; \cal{H} in \cal{X})

#' @title Size of human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H <- function(t, y, pars) {
  UseMethod("F_H", pars$Hpar)
}

#' @title Derivatives of demographic changes in human populations
#' @description This method dispatches on the type of `y`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return see help pages for specific methods
#' @export
dHdt <- function(t, y, pars){
  UseMethod("dHdt", y)
}

#' @title A function that computes the birth rate for human populations
#' @description This method dispatches on the type of `y`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return see help pages for specific methods
#' @export
Births <- function(t, y, pars){
  UseMethod("Births", y)
}

#' @title Add indices for human population denominators to parameter list
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars a [list]
#' @return none
#' @export
make_indices_H <- function(pars) {
  UseMethod("make_indices_H", pars$Hpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_H <- function(pars) {
  UseMethod("get_inits_H", pars$Hpar)
}

#' @title Make parameters for null human demography model
#' @param pars a [list]
#' @param H size of human population in each strata
#' @param residence is a vector describing patch residency
#' @param searchWts is a vector describing blood feeding search weights
#' @param TaR is a matrix describing time spent among patches
#' @return none
#' @export
make_parameters_demography_null <- function(pars, H, residence, searchWts, TaR) {
  stopifnot(length(H) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c("static")
  Hpar$H <- H

  class(Hpar$H) <- "static"
  Hpar$residence <- residence
  Hpar$wts_f <- searchWts
  Hpar$TaR <- TaR

  birthF <- "null"
  class(birthF) <- "null"
  Hpar$birthF <- birthF
  Hpar$Hmatrix <- diag(0, length(H))

  pars$Hpar <- Hpar
  pars$nStrata <- length(H)

  return(pars)
}


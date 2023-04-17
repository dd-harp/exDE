# generic methods for demography (nested within human; \cal{H} in \cal{X})

#' @title Size of human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H <- function(t, y, pars) {
  UseMethod("F_H", pars$Hpar)
}

#' @title Derivatives of demographic changes in human populations
#' @description This method dispatches on the type of `y`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return see help pages for specific methods
#' @export
dHdt <- function(t, y, pars){
  UseMethod("dHdt", y)
}

#' @title A function that computes the birth rate for human populations
#' @description This method dispatches on the type of `y`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return see help pages for specific methods
#' @export
Births <- function(t, y, pars){
  UseMethod("Births", y)
}

#' @title Add indices for human population denominators to parameter list
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @return none
#' @export
make_indices_H <- function(pars) {
  UseMethod("make_indices_H", pars$Hpar)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_H <- function(pars) {
  UseMethod("get_inits_H", pars$Hpar)
}

#' @title Make parameters for null human demography model
#' @param pars an [environment]
#' @param H size of human population in each strata
#' @param membershipH is a vector describing patch residency
#' @param searchWtsH is a vector describing blood feeding search weights
#' @param TimeSpent is a matrix describing time spent among patches
#' @return none
#' @export
make_parameters_demography_null <- function(pars, H, membershipH, searchWtsH, TimeSpent) {
  stopifnot(length(H) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c("static")
  Hpar$H <- H
  class(Hpar$H) <- "static"
  Hpar$membershipH <- membershipH
  Hpar$searchWtsH <- searchWtsH
  Hpar$wf <- searchWtsH
  Hpar$TimeSpent <- TimeSpent
  Hpar$Psi <- TimeSpent
  birthF <- "null"
  class(birthF) <- "null"
  Hpar$birthF <- birthF
  Hpar$Hmatrix <- diag(0, length(H))
  pars$Hpar <- Hpar
  pars$nStrata <- length(H)
  pars$beta <- compute_beta(H, searchWtsH, TimeSpent)
  pars$beta_lag <- pars$beta
  return(pars)
}


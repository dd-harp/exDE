# null model of \mathcal{H}; constant for all time

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the null model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length 0
#' @export
Births.null <- function(t, y, pars){
  if(class(y)=="null") numeric(0) else 0*y
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the null model.
#' @inheritParams dHdt
#' @return a [numeric] vector of 0s or of length 0
#' @export
dHdt.null <- function(t, y, pars){
  if(class(y)=="null") numeric(0) else 0*y
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
  class(Hpar) <- c('static')
  Hpar$H <- H
  class(Hpar$H) <- "null"
  Hpar$membershipH <- membershipH
  Hpar$searchWtsH <- searchWtsH
  Hpar$TimeSpent <- TimeSpent
  birthF <- "null"
  class(birthF) <- "null"
  Hpar$birthF <- birthF
  Hmatrix <- "null"
  class(Hmatrix) <- "null"
  Hpar$Hmatrix <- Hmatrix
  pars$Hpar <- Hpar
  pars$nStrata <- length(H)
  pars$beta <- compute_beta(H, searchWtsH, TimeSpent)
  pars$beta_lag <- pars$beta
  return(pars)
}

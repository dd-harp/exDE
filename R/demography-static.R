# static model of \mathcal{H}; constant for all time

#' @title Size of human population denominators
#' @description Implements [F_H] for the static model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.static <- function(t, y, pars) {
  pars$Hpar$H
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [Births] when `y` is static
#' @inheritParams Births
#' @return a [numeric] vector of length 0
#' @export
Births.static <- function(t, y, pars){
  numeric(0)
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] when `y` is static
#' @inheritParams dHdt
#' @return a [numeric] vector of length 0
#' @export
dHdt.static <- function(t, y, pars){
  numeric(0)
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_indices_H] for static model.
#' @inheritParams make_indices_H
#' @return none
#' @export
make_indices_H.static <- function(pars) {
  pars$H_ix <- integer(0)
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_H.static<- function(pars){
  return(numeric(0))
}

#' @title Make parameters for static human demography model
#' @param pars an [environment]
#' @param H size of human population in each strata
#' @param membershipH is a vector describing patch residency
#' @param searchWtsH is a vector describing blood feeding search weights
#' @param TimeSpent is a matrix describing time spent among patches
#' @param birthFpars setup to dispatch and compute `F_birth`
#' @param Hmatrix does a set of state transitions
#' @param birthsXstrata distributes births to the youngest strata
#' @return none
#' @export
make_parameters_demography_static <- function(pars, H, membershipH, searchWtsH, TimeSpent,
                                              birthFpars, Hmatrix, birthsXstrata) {
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
  Hpar$birthFpars <- birthFpars
  Hpar$birthXstrata <- birthsXstrata
  Hpar$Hmatrix <- Hmatrix
  pars$Hpar <- Hpar
  pars$nStrata <- length(H)
  pars$beta <- compute_beta(H, searchWtsH, TimeSpent)
  pars$beta_lag <- pars$beta
  return(pars)
}

# basic driven model of demography

#' @title Size of human population denominators
#' @description Implements [F_H] for the forced (basic) model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.basic <- function(t, y, pars) {
  y[pars$H_ix]
}

#' @title Size of lagged human population denominators
#' @description Implements [F_H_lag] for forced (basic) model.
#' @inheritParams F_H_lag
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H_lag.basic <- function(t, y, pars, lag) {
  pars$Hpar$H
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the forced (basic) model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata` or of length 0
#' @export
dHdt.basic <- function(t, y, pars){
  pars$Hpar$deaths%*%y
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the forced (basic) model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata` or of length 0
#' @export
Births.basic <- function(t, y, pars){with(pars$Hpar,{
  F_birth(t, y, pars)*birthsXstrata
})}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_indices_H] for forced (basic) model.
#' @inheritParams make_indices_H
#' @return none
#' @export
make_indices_H.basic <- function(pars) {
  pars$H_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$H_ix, 1)
  return(pars)
}

#' @title Make parameters for forced (basic) human demography model
#' @param pars an [environment]
#' @param H a function taking a single argument `t` and returning a vector of length
#' `nStrata`.
#' @return none
#' @export
make_parameters_demography_basic <- function(pars, H, F_birth, birthrate, deathrate, birthsXstrata) {
  stopifnot(length(H) == pars$nStrata)
  stopifnot(length(deathrate) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c('basic')
  Hpar$H <- H
  Hpar$F_birth <- F_birth
  Hpar$birthrate <- birthrate
  Hpar$birthXstrata <- birthsXstrata
  Hpar$deaths <- diag(-deathrate)
  pars$Hpar <- Hpar
  return(pars)
}

#' @title Attach initial values to parameter list
#' @description Implements [attach_inits_H] for the basic model.
#' @inheritParams attach_inits_H
#' @return none
#' @export
attach_inits_H.basic <- function(pars) {
  pars$Hpar$H0 = pars$Hpar$inits$H0
  return(pars)
}

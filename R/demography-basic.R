# basic model of demography with births and deaths

#' @title Size of human population denominators
#' @description Implements [F_H] for the basic demography model
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.basic <- function(t, y, pars) {
  y[pars$H_ix]
}

#' @title Size of lagged human population denominators
#' @description Implements [F_H_lag] for the basic demography model
#' @inheritParams F_H_lag
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H_lag.basic <- function(t, y, pars, lag) {
  pars$Hpar$H
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the basic demography model
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata` or of length 0
#' @export
dHdt.basic <- function(t, y, pars){
  pars$Hpar$deaths%*%y
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the basic demography model
#' @inheritParams Births
#' @return a [numeric] vector of length `nStrata` or of length 0
#' @export
Births.basic <- function(t, y, pars){with(pars$Hpar,{
  F_birth(t, y, pars)*birthsXstrata
})}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_indices_H] for the basic demography model.
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
#' @param membershipH is a vector describing patch residency
#' @param searchWtsH is a vector describing blood feeding search weights
#' @param TimeSpent is a matrix describing time spent among patches
#' @param F_birth is function that returns the birth rate as a function of time
#' @param birthrate is a parameter describing the current birthrate for each stratum
#' @param deathrate is a parameter describing the current deathrate for each stratum
#' @param birthsXstrata is a parameter describing which strata get births
#' @return none
#' @export
make_parameters_demography_basic <- function(pars, H, membershipH, searchWtsH, TimeSpent,
                                             F_birth, birthrate, deathrate, birthsXstrata) {
  stopifnot(length(H) == pars$nStrata)
  stopifnot(length(deathrate) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c('basic')
  Hpar$H <- H
  Hpar$membershipH <- membershipH
  Hpar$searchWtsH <- searchWtsH
  Hpar$TimeSpent <- TimeSpent
  Hpar$F_birth <- F_birth
  Hpar$birthrate <- birthrate
  Hpar$birthXstrata <- birthsXstrata
  Hpar$deaths <- diag(-deathrate)
  pars$Hpar <- Hpar
  return(pars)
}


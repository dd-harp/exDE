# dynamic model of demography with births and deaths

#' @title Size of human population denominators
#' @description Implements [F_H] for the dynamic demography model
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.dynamic <- function(t, y, pars) {
  y[pars$H_ix]
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_indices_H] for the dynamic demography model.
#' @inheritParams make_indices_H
#' @return none
#' @export
make_indices_H.dynamic <- function(pars) {
  pars$H_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$H_ix, 1)
  return(pars)
}

#' @title Parse the output of deSolve and return H for models with dynamic demography
#' @description Implements [parse_deout_H] for models with dynamic demography
#' @inheritParams parse_deout_H
#' @return none
#' @export
parse_deout_H.dynamic <- function(varslist, deout, pars) {
  varslist$H = deout[,pars$Hpar$H_ix+1]
  return(varslist)
}

#' @title Update inits for the static human demography model
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_H.dynamic <- function(pars, y0) {
  pars$Hpar$H = y0[pars$Hpar$H_ix]
  return(pars)
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_H.dynamic <- function(pars){
  return(pars$Hpar$H)
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [Births] when `y` is numeric
#' @inheritParams Births
#' @return a [numeric] vector of length `nStrata`
#' @export
Births.numeric <- function(t, y, pars){
  F_births(t, y, pars)*pars$Hpar$birthsXstrata
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] when `y` is numeric
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata`
#' @export
dHdt.numeric <- function(t, y, pars){
  pars$Hpar$Hmatrix %*% y
}

#' @title Make parameters for forced (dynamic) human demography model
#' @param pars a [list]
#' @param H a function taking a single argument `t` and returning a vector of length
#' `nStrata`.
#' @param residence is a vector describing patch residency
#' @param searchWts is a vector describing blood feeding search weights
#' @param TaR is a matrix describing time spent among patches
#' @param birthFpars setup to dispatch and compute `F_birth`
#' @param Hmatrix does a set of state transitions
#' @param birthsXstrata distributes births to the youngest strata
#' @return none
#' @export
make_parameters_demography_dynamic <- function(pars, H, residence, searchWts, TaR,
                                             birthFpars, Hmatrix, birthsXstrata) {
  stopifnot(length(H) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c('dynamic')
  Hpar$H <- H
  Hpar$residence <- residence
  Hpar$wts_f <- searchWts
  Hpar$TaR <- TaR
  Hpar$birthFpars <- birthFpars
  Hpar$birthXstrata <- birthsXstrata
  Hpar$Hmatrix <- Hmatrix
  pars$Hpar <- Hpar
  pars$nStrata <- length(H)
  return(pars)
}



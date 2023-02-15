# null model of \mathcal{H}; constant for all time

#' @title Size of human population denominators
#' @description Implements [F_H] for the null model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.null <- function(t, y, pars) {
  pars$Hpar$H
}

#' @title Size of lagged human population denominators
#' @description Implements [F_H_lag] for the null model.
#' @inheritParams F_H_lag
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H_lag.null <- function(t, y, pars, lag) {
  pars$Hpar$H
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the null model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length 0
#' @export
Births.null <- function(t, y, pars){
  if(class(y) == 'param') numeric(0) else 0*y
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for the null model.
#' @inheritParams dHdt
#' @return a [numeric] vector of length 0
#' @export
dHdt.null <- function(t, y, pars){
  if(class(y) == 'param') numeric(0) else 0*y
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_indices_H] for null model.
#' @inheritParams make_indices_H
#' @return none
#' @export
make_indices_H.null <- function(pars) {
  pars$H_ix <- integer(0)
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars an [environment]
#' @return none
#' @export
get_inits_H.null<- function(pars){
  return(numeric(0))
}

#' @title Make parameters for null human demography model
#' @param pars an [environment]
#' @param H size of human population in each strata
#' @return none
#' @export
make_parameters_demography_null <- function(pars, H, membershipH, searchWtsH, TimeSpent) {
  stopifnot(length(H) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c('null')
  Hpar$H <- H
  class(Hpar$H) <- 'param'
  Hpar$membershipH <- membershipH
  Hpar$searchWtsH <- searchWtsH
  Hpar$TimeSpent <- TimeSpent
  class(Hpar$TimeSpent) <- 'static'
  pars$Hpar <- Hpar
  pars$nStrata = length(H)
  return(pars)
}

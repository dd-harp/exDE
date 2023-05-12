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
#' @param pars a [list]
#' @return none
#' @export
get_inits_H.static<- function(pars){
  return(numeric(0))
}

#' @title Make parameters for static human demography model
#' @param pars a [list]
#' @param H size of human population in each strata
#' @param residence is a vector describing patch residency
#' @param searchWts is a vector describing blood feeding search weights
#' @param TaR is a matrix describing time spent among patches
#' @param birthFpars setup to dispatch and compute `F_birth`
#' @param Hmatrix does a set of state transitions
#' @param birthsXstrata distributes births to the youngest strata
#' @return none
#' @export
make_parameters_demography_static <- function(pars, H, residence, searchWts, TaR,
                                              birthFpars, Hmatrix, birthsXstrata) {
  stopifnot(length(H) == pars$nStrata)
  Hpar <- list()
  class(Hpar) <- c("static")
  Hpar$H <- H
  class(Hpar$H) <- "static"
  Hpar$residence <- residence
  Hpar$wts_f <- searchWts
  Hpar$TaR <- TaR

  Hpar$birthFpars <- birthFpars
  Hpar$birthXstrata <- birthsXstrata
  Hpar$Hmatrix <- Hmatrix
  pars$Hpar <- Hpar
  pars$nStrata <- length(H)
  pars$beta <- compute_beta(H, searchWts, TaR)
  pars$beta_lag <- pars$beta
  return(pars)
}

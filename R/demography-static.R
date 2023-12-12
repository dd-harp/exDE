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


#' @title Update inits for the static human demography model
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_H.static <- function(pars, y0) {
  return(pars)
}

#' @title Parse the output of deSolve and return variables for models where H is a parameter
#' @description Implements [parse_deout_H] for models with constant denominators
#' @inheritParams parse_deout_H
#' @return a [list]
#' @export
parse_deout_H.static <- function(deout, pars) {
  H = pars$Hpar$H
  H = matrix(H, nrow = length(deout[,1]), ncol = length(H), byrow=T)
  return(list(H=H))
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
  Hpar$residence <- checkIt(residence, pars$nStrata, F)
  Hpar$wts_f <- checkIt(searchWts, pars$nStrata, F)
  Hpar$rbr <- searchWts*sum(H)/sum(searchWts*H)
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

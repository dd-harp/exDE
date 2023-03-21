# specialized methods for the aquatic mosquito trace model

#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the trace (forced emergence) model.
#' @inheritParams F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.trace <- function(t, y, pars) {
  with(pars$Lpar, Lambda*Lt(t, pars))
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the trace (forced emergence) model.
#' @inheritParams dLdt
#' @return a [numeric] vector
#' @export
dLdt.trace <- function(t, y, pars, eta) {
  numeric(0)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_L] for trace (forced emergence) model.
#' @inheritParams make_indices_L
#' @return none
#' @export
make_indices_L.trace <- function(pars) {
  pars$L_ix <- integer(0)
  return(pars)
}

#' @title Make parameters for trace aquatic mosquito model
#' @param pars an [environment]
#' @param Lambda vector of mean emergence rates from each aquatic habitat
#' @param Lt is a [function] of the form Lt(t,pars) that computes temporal fluctuations
#' @return none
#' @export
make_parameters_L_trace <- function(pars, Lambda, Lt="one") {
  stopifnot(is.numeric(Lambda))
  Lpar <- list()
  class(Lpar) <- 'trace'
  Lpar$Lambda <- Lambda
  if(Lt == "one"){
    Lpar$Lt = function(t, pars){1}
  }
  pars$Lpar <- Lpar
  return(pars)
}

#' @title Make parameters for trace aquatic mosquito model
#' @param pars an [environment]
#' @param L0 is set to NULL for the trace model
#' @return none
#' @export
make_inits_L_trace<- function(pars, L0=NULL) {
  pars$Linits = numeric(0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_L] for the GeRM model.
#' @inheritParams get_inits_L
#' @return none
#' @export
get_inits_L.trace <- function(pars){
  numeric(0)
}


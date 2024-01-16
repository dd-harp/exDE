#' @title Exposure and Infection
#' @description This function translates the local daily entomological
#' inoculation rate (dEIR), computed as one of the transmission terms
#' into a measure of the daily force of infection (dFoI). The daily FoI
#' is the sum of two terms: 1) a function [F_foi] computes the local dFoI;
#' 2) a function [travel_malaria] that computes the FoI resulting from
#' exposure while traveling.
#' @param t the time
#' @param y the state variables
#' @param pars the model object as a [list]
#' @return the function modifies **pars** and returns it: the computed FoI are stored as `pars$FoI`
#' @export
Exposure <- function(t, y, pars){

  for(i in 1:pars$nHosts){
    b = F_b(y, pars, i)
    pars$FoI[[i]] = F_foi(pars$EIR[[i]], b, pars) + travel_malaria(t, pars)
  }
  return(pars)
}

#' @title A model for daily FoI as a function of the daily EIR.
#' @description This function compute the daily local FoI as a function
#' of the daily EIR and effects of partial immunity.
#' It is computed based on a probabilistic model of exposure and
#' possibly including environmental heterogeneity. If a model of human / host
#' infection has terms that describe partial immunity, *e.g.* affecting
#' pre-erythrocytic immunity to malaria called by [F_b], those effects are implemented here.
#' The method dispatches on the type of `pars$FOIpar`.
#' @param eir the daily eir for each stratum
#' @param b the probability of infection, per bite
#' @param pars the model object as a [list]
#' @return the daily, local force of infection as a [numeric] vector
#' @export
F_foi <- function(eir, b, pars){
  UseMethod("F_foi", pars$FOIpar)
}

#' @title A Poisson model for the daily local FoI as a function of the daily EIR.
#' @description Implements [F_foi] for a Poisson model
#' @inheritParams F_foi
#' @return a [numeric] vector of length `nStrata`
#' @export
F_foi.pois <- function(eir, b, pars){
  b*eir
}

#' @title The daily FoI as a function of the daily EIR under a negative binomial model of exposure.
#' @description Implements [F_foi] for a negative binomial model
#' @inheritParams F_foi
#' @return a [numeric] vector of length `nStrata`
#' @export
F_foi.nb <- function(eir, b, pars){
  log(1 + b*eir/pars$FOIpar$sz)*pars$FOIpar$sz
}

#' @title Make parameters for the null model of exposure
#' @param pars a [list]
#' @return none
#' @export
setup_exposure_pois <- function(pars) {
  FOIpar <- list()
  class(FOIpar) <- 'pois'
  pars$FOIpar <- FOIpar
  return(pars)
}

#' @title Make parameters for the null model of exposure
#' @param pars a [list]
#' @param sz the size parameter, as in dnbinom(mu=mu, size=size)
#' @return none
#' @export
setup_exposure_nb <- function(pars, sz) {
  FOIpar <- list()
  class(FOIpar) <- 'nb'
  FOIpar$sz = sz
  pars$FOIpar <- FOIpar
  return(pars)
}

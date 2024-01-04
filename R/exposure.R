#' @title A model for exposure. The function `F_b` must be define
#' @param t the time
#' @param y the variables
#' @param pars a [list]
#' @return [list]
#' @export
Exposure <- function(t, y, pars){
  for(i in 1:pars$nHosts){
    b = F_b(y, pars, i)
    pars$FoI[[i]] = F_foi(pars$EIR[[i]], b, pars) + travel_malaria(t, pars)
  }
  return(pars)
}

#' @title A model for exposure
#' @description This method dispatches on the type of `pars$FOIpar`.
#' @param eir the daily eir
#' @param b the probability of infection, per bite
#' @param pars a [list]
#' @return see help pages for specific methods
#' @export
F_foi <- function(eir, b, pars){
  UseMethod("F_foi", pars$FOIpar)
}

#' @title Null human population births
#' @description Implements [F_foi] for a Poisson model
#' @inheritParams F_foi
#' @return a [numeric] vector of length `nStrata`
#' @export
F_foi.pois <- function(eir, b, pars){
  b*eir
}

#' @title Null human population births
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

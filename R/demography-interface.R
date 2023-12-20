# generic methods for demography (nested within human; \cal{H} in \cal{X})

#' @title A utility to set up Hpar
#' @param pars a [list]
#' @param HPop a [numeric] vector of population densities
#' @param residence is the patch where each stratum resides
#' @param searchWts is the blood feeding search weight for each stratum
#' @param Hopts a [list] to overwrite default values
#' @return a [list]
#' @export
setup_Hpar = function(pars, HPop=1000, residence=1, searchWts=1, Hopts=NULL){
  with(Hopts,{
    with(pars, if(!exists("nStrata")) pars$nStrata = length(HPop))

    Hpar = list()
    Hpar$H = HPop

    Hpar$residence <- checkIt(residence, pars$nStrata, F)
    Hpar$wts_f <- checkIt(searchWts, pars$nStrata, F)
    Hpar$rbr <- searchWts*sum(HPop)/sum(searchWts*HPop)

    Births <- "zero"
    class(Births) <- "zero"
    Hpar$Births <- Births

    dH <- "zero"
    class(dH) <- "zero"
    Hpar$dH <- dH

    Hpar$residence = residence
    Hpar$wts_f = searchWts

    pars$Hpar <- Hpar
    return(pars)
})}


#' @title Derivatives of demographic changes in human populations
#' @description This method dispatches on the type of pars$Hpar$dH
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return see help pages for specific methods
#' @export
dHdt <- function(t, y, pars){
  UseMethod("dHdt", pars$Hpar$dH)
}

#' @title A function that computes the birth rate for human populations
#' @description This method dispatches on the type of pars$Hpar$Births
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return see help pages for specific methods
#' @export
Births <- function(t, y, pars){
  UseMethod("Births", pars$Hpar$Births)
}

#' @title Make parameters for null human demography model
#' @param pars a [list]
#' @param H size of human population in each strata
#' @param residence is a vector describing patch residency
#' @param searchWts is a vector describing blood feeding search weights
#' @param TaR is a matrix describing time spent among patches
#' @return none
#' @export
make_parameters_demography_null <- function(pars, H, residence, searchWts, TaR) {
  stopifnot(length(H) == pars$nStrata)
  Hpar <- list()
  Hpar$H <- H

  Hpar$residence <- checkIt(residence, pars$nStrata, F)
  Hpar$wts_f <- checkIt(searchWts, pars$nStrata, F)
  Hpar$rbr <- searchWts*sum(H)/sum(searchWts*H)
  Hpar$TaR <- TaR

  Births <- "zero"
  class(Births) <- "zero"
  Hpar$Births <- Births

  dH <- "zero"
  class(dH) <- "zero"
  Hpar$dH <- dH

  pars$Hpar <- Hpar
  pars$nStrata <- length(H)

  return(pars)
}



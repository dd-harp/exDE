# generic methods for demography (nested within human; \cal{H} in \cal{X})

#' @title A utility to set up Hpar
#' @param pars a [list]
#' @param i the host species index
#' @param HPop a [numeric] vector of population densities
#' @return a [list]
#' @export
setup_Hpar_static = function(pars, i, HPop=1000){

  Hpar = list()
  class(Hpar) <- "static"
  Hpar$H = HPop
  Hpar$nStrata = length(HPop)

  Bf <- "zero"
  class(Bf) <- "zero"
  Hpar$Bf <- Bf

  dH <- "zero"
  class(dH) <- "zero"
  Hpar$dH <- dH

  pars$Hpar[[i]] <- Hpar

  return(pars)
}


#' @title Derivatives of demographic changes in human populations
#' @description This method dispatches on the type of pars$Hpar$dH
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param i the index of the host species
#' @return see help pages for specific methods
#' @export
dHdt <- function(t, y, pars, i){
  UseMethod("dHdt", pars$Hpar[[i]]$dH)
}

#' @title A function that computes the birth rate for human populations
#' @description This method dispatches on the type of pars$Hpar$Births
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param i the index of the host species
#' @return see help pages for specific methods
#' @export
Births <- function(t, y, pars, i){
  UseMethod("Births", pars$Hpar[[i]]$Bf)
}

#' @title Make parameters for null human demography model
#' @param pars a [list]
#' @param H size of human population in each strata
#' @return none
#' @export
make_parameters_demography_null <- function(pars, H) {
  stopifnot(length(H) == pars$nStrata)
  Hpar <- list()
  Hpar$H <- H
  Hpar$nStrata <- length(H)

  Bf <- "zero"
  class(Bf) <- "zero"
  Hpar$Bf <- Bf

  dH <- "zero"
  class(dH) <- "zero"
  Hpar$dH <- dH

  pars$Hpar[[1]] <- Hpar

  return(pars)
}



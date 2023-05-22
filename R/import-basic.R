# specialized methods for the basic model of parasite / pathogen importation

#' @title Visitors, a basic model
#' @description Implements [Visitors] for the basic model of importation (do nothing)
#' @inheritParams Visitors
#' @return a named [list]
#' @export
Visitors.basic <- function(t, pars) {
  pars$x_visitors =  with(pars$Ipar, xv*xt(t, pars))
  pars$Visitors =  with(pars$Ipar, Vm*Vt(t, pars))
  return(pars)
}

#' @title Parasite / pathogen importation, the basic model
#' @description Implements [kappa_local] for the basic model (do nothing)
#' @inheritParams kappa_local
#' @return kappa, a numeric [vector]
#' @export
kappa_local.basic<- function(kappa, pars) {
  with(pars, local_frac*kappa + (1-local_frac)*x_visitors)
  return(kappa)
}

#' @title Parasite / pathogen importation, the basic model
#' @description Implements [fqZ_local] for the basic model (do nothing)
#' @inheritParams fqZ_local
#' @return fqZ, a numeric [vector]
#' @export
fqZ_local.basic<- function(fqZ, pars) {
  return(pars$local_frac*fqZ)
}

#' @title Make parameters for the basic model for visitors
#' @param pars a [list]
#' @param local_frac the fraction of bites taken on residents
#' @param xv the mean probability a bite on a visitor infects a mosquito
#' @param xt a function to force x_visitors over time
#' @param Vm is the mean availability of visitors for blood feeding
#' @param Vt a function to force availability of Visitors over time
#' @return none
#' @export
setup_visitors_basic <- function(pars, local_frac=1, xv=0, xt=NULL, Vm=0, Vt = NULL) {
  Ipar <- list()
  class(Ipar) <- 'basic'
  pars$Ipar <- Ipar
  pars$Ipar$xv = xv
  if(is.null(xt)) xt = function(t, pars){1}
  pars$local_frac = 0
  pars$Ipar$Vm = Vm
  if(is.null(Vt)) Vt = function(t, pars){1}
  return(pars)
}


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

#' @title A function to set up malaria importation
#' @description Implements [setup_MalariaImportation] for the basic model
#' @inheritParams setup_MalariaImportation
#' @return a [list]
#' @export
setup_MalariaImportation.null = function(pars, IMname, IMopts =list()){
  Ipar <- list()
  class(Ipar) <- 'basic'
  pars$Ipar <- Ipar
  pars = setup_Visitors_basic(pars, IMopts)
  pars = setup_Travel_basic(pars, IMopts)
  return(pars)
}

#' @title Make parameters for the basic model for visitors
#' @param pars a [list]
#' @param IMopts a [list]
#' @param local_frac the fraction of bites taken on residents
#' @param xv the mean probability a bite on a visitor infects a mosquito
#' @param xt a function to force x_visitors over time
#' @param Vm is the mean availability of visitors for blood feeding
#' @param Vt a function to force availability of Visitors over time
#' @return none
#' @export
setup_Visitors_basic <- function(pars, IMopts, local_frac=1, xv=0, xt=NULL, Vm=0, Vt = NULL) {with(IMopts,{

  pars$local_frac = 1

  pars$Ipar$xv = xv
  if(is.null(xt)) xt = function(t, pars){1}
  pars$Ipar$xt = xt

  pars$Ipar$Vm = Vm
  if(is.null(Vt)) Vt = function(t, pars){1}
  pars$Ipar$Vt = Vt

  return(pars)
})}

#' @title Make parameters for the basic model for Travel Malaria
#' @param pars a [list]
#' @param TravelOpts a list, overwrites default values
#' @param delta the mean FoI for travel malaria
#' @param delta_t the temporal signal (should have a mean of 1)
#' @return none
#' @export
setup_Travel_basic <- function(pars, TravelOpts = list(), delta=0, delta_t=NULL) {with(TravelOpts,{
  pars$Ipar$delta = delta
  if(is.null(delta_t)) delta_t = function(t, pars){1}
  pars$Ipar$delta_t = delta_t
  return(pars)
})}

#' @title travel_foi, a null model
#' @description Implements [travel_foi] for the null model (do nothing)
#' @inheritParams travel_foi
#' @return a [numeric]
#' @export
travel_foi.basic <- function(t, pars) {
  reteurn(with(pars$Ipar, delta*delta_t(t, pars)))
}



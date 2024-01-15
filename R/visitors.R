# generic methods for parasite / pathogen importation by visitors

#' @title Visitors
#' @description This method dispatches on the type of `pars$VISITORS`.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
Visitors <- function(t, pars) {
  UseMethod("Visitors", pars$VISITORS)
}

#' @title Visitors, a static model
#' @description Implements [Visitors] for the static model (do nothing)
#' @inheritParams Visitors
#' @return a [list]
#' @export
Visitors.static <- function(t, pars) {
  return(pars)
}

#' @title Make parameters for the static model visitors (no visitors)
#' @param pars a [list]
#' @param local_frac is the fraction of humans / hosts that are not visitors
#' @param Visitors is the availability of visitors
#' @param x_visitors is the net infectiousness of the visitors
#' @return [list]
#' @export
setup_visitors_static <- function(pars, local_frac=1, Visitors=0, x_visitors=0) {

  VISITORS <- list()
  class(VISITORS) <- "static"
  pars$VISITORS <- VISITORS

  pars$vars$local_frac = list()
  pars$vars$local_frac[[1]] = local_frac
  pars$vars$Visitors = list()
  pars$vars$Visitors[[1]]  = Visitors
  pars$vars$x_visitors = list()
  pars$vars$x_visitors[[1]]  = x_visitors

  return(pars)
}


#' @title Visitors, the basic model
#' @description Implements [Visitors] for the basic model for Visitors
#' @inheritParams Visitors
#' @return a [list]
#' @export
Visitors.basic <- function(t, pars) {
  pars$vars$x_visitors =  with(pars$VISITORS, x_scale*xt(t, pars))
  pars$vars$Visitors =  with(pars$VISITORS, V_scale*Vt(t, pars))
  return(pars)
}

#' @title Make parameters and functions for the basic model for visitors
#' @param pars a [list]
#' @param IMopts a [list]
#' @param x_scale a non-negative numeric value to set the mean for x_visitors
#' @param xt a function to change the pattern for x_visitors over time
#' @param V_scale a non-negative numeric value to set the mean availability of Visitors
#' @param Vt a function to set the temporal pattern for availability of Visitors
#' @return [list]
#' @export
setup_visitors_basic <- function(pars, IMopts, x_scale = 0, xt = NULL, V_scale = 0, Vt = NULL) {with(IMopts,{

  pars$VISITORS$x_scale = x_scale
  if(is.null(xt)) xt = function(t, pars){1}
  pars$VISITORS$xt = xt

  pars$VISITORS$V_scale = V_scale
  if(is.null(Vt)) Vt = function(t, pars){1}
  pars$VISITORS$Vt = Vt

  return(pars)
})}

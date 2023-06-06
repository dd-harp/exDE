# specialized methods for a basic visitors model

#' @title Visitors, the basic model
#' @description Implements [Visitors] for the basic model for Visitors
#' @inheritParams Visitors
#' @return a named [list]
#' @export
Visitors.basic <- function(t, pars) {
  pars$x_visitors =  with(pars$Ipar, x_scale*xt(t, pars))
  pars$Visitors =  with(pars$Ipar, V_scale*Vt(t, pars))
  return(pars)
}

#' @title Make parameters and functions for the basic model for visitors
#' @param pars a [list]
#' @param IMopts a [list]
#' @param x_scale a non-negative numeric value to set the mean for x_visitors
#' @param xt a function to change the pattern for x_visitors over time
#' @param V_scale a non-negative numeric value to set the mean availability of Visitors
#' @param Vt a function to set the temporal pattern for availability of Visitors
#' @return none
#' @export
setup_visitors_basic <- function(pars, IMopts, x_scale = 0, xt = NULL, V_scale = 0, Vt = NULL) {with(IMopts,{

  pars$Ipar$x_scale = x_scale
  if(is.null(xt)) xt = function(t, pars){1}
  pars$Ipar$xt = xt

  pars$Ipar$V_scale = V_scale
  if(is.null(Vt)) Vt = function(t, pars){1}
  pars$Ipar$Vt = Vt

  return(pars)
})}

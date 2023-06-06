
#' @title Set the availability of other
#' @description This method dispatches on the type of `pars$OBHpar`. It should
#' set the values of other availability at time `t`
#' @param t current simulation time
#' @param pars a [list]
#' @return a [list]
#' @export
F_otherblood <- function(t, pars) {
  UseMethod("F_otherblood", pars$OBHpar)
}

#' @title Null model for other blood hosts
#' @description Implements [F_otherblood] for a null model
#' @inheritParams F_otherblood
#' @return a [list]
#' @export
F_otherblood.null <- function(t, pars){
  pars
}

#' @title Static model for other blood hosts
#' @description Implements [F_otherblood] for a static model
#' @inheritParams F_otherblood
#' @return a [list]
#' @export
F_otherblood.static <- function(t, pars){
  pars
}

#' @title A basic forcing model for other blood hosts
#' @description Implements [F_otherblood] for a forced model
#' @inheritParams F_otherblood
#' @return a [list]
#' @export
F_otherblood.forced <- function(t, pars){
  with(pars$OBHpar,{
    return(otherblood_scale*otherblood_t(t, pars))
  })
}

#' @title Set up forcing for the availability of other blood hosts
#' @param pars a [list]
#' @param otherblood_Opts a list, overwrites default values
#' @param otherblood_scale a non-negative numeric value to scale the mean availability of other blood hosts
#' @param otherblood_t a function that specifies the temporal pattern for the availability of other blood hosts
#' @return none
#' @export
setup_otherblood_forced <- function(pars, otherblood_Opts = list(), otherblood_scale=0, otherblood_t=NULL) {with(otherblood_Opts,{

  OBHpar <- list()
  class(OBHpar) <- 'foi'
  pars$OBHpar <- OBHpar

  pars$OBHpar$otherblood_scale = otherblood_scale
  if(is.null(otherblood_t)) otherblood_t = function(t, pars){1}
  pars$OBHpar$otherblood_t = otherblood_t

  return(pars)
})}



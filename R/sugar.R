
#' @title Set the availability of sugar
#' @description This method dispatches on the type of `pars$SGRpar`. It should
#' set the values of sugar availability at time `t`
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sugar <- function(t, pars) {
  UseMethod("F_sugar", pars$SGRpar)
}

#' @title Null model for sugar
#' @description Implements [F_sugar] for a null model
#' @inheritParams F_sugar
#' @return a [numeric] vector of length `nStrata`
#' @export
F_sugar.null <- function(t, pars){
  pars
}

#' @title Basic model for sugar
#' @description Implements [F_sugar] for a basic model
#' @inheritParams F_sugar
#' @return a [numeric] vector of length `nStrata`
#' @export
F_sugar.static <- function(t, pars){
  pars
}

#' @title A basic forcing model for sugar
#' @description Implements [F_sugar] for a forced model
#' @inheritParams F_sugar
#' @return a [list]
#' @export
F_sugar.forced <- function(t, pars){
  with(pars$SGRpar,{
    return(sugar_scale*sugar_t(t, pars))
  })
}

#' @title Set up forcing for the availability of sugar
#' @param pars a [list]
#' @param sugar_Opts a list, overwrites default values
#' @param sugar_scale a non-negative numeric value to scale the mean availability of sugar
#' @param sugar_t a function that specifies the temporal pattern for the availability of sugar
#' @return none
#' @export
setup_sugar_forced <- function(pars, sugar_Opts = list(), sugar_scale=0, sugar_t=NULL) {with(sugar_Opts,{

  SGRpar <- list()
  class(SGRpar) <- 'foi'
  pars$SGRpar <- SGRpar

  pars$SGRpar$sugar_scale = sugar_scale
  if(is.null(sugar_t)) sugar_t = function(t, pars){1}
  pars$SGRpar$sugar_t = sugar_t

  return(pars)
})}

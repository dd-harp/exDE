#' @title Set the values for habitat
#' @description This method dispatches on the type of `pars$QHABpar`. It should
#' set the values of habitat at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
F_habitats <- function(t, y, pars) {
  UseMethod("F_habitat", pars$QHABpar)
}

#' @title A null model for habitat
#' @description Implements [F_habitats] for a static model
#' @inheritParams F_habitats
#' @return a [list]
#' @export
F_habitats.null <- function(t, y, pars){
  pars
}

#' @title Set up the null model for habitat dynamics
#' @param pars a [list]
#' @return none
#' @export
setup_habitats_null<- function(pars){
  QHABpar <- list()
  class(QHABpar) <- 'null'
  pars$QHABpar <- QHABpar
  return(pars)
}

#' @title Static model for habitat
#' @description Implements [F_habitats] for a static model
#' @inheritParams F_habitats
#' @return a [list]
#' @export
F_habitats.static <- function(t, y, pars){
  pars
}

#' @title Set up the static model for habitat dynamics
#' @param pars a [list]
#' @return none
#' @export
setup_habitats_static<- function(pars){
  QHABpar <- list()
  class(QHABpar) <- 'static'
  pars$QHABpar <- QHABpar
  return(pars)
}


#' @title A basic forcing model for habitat availability
#' @description Implements [F_habitats] for a forced model
#' @inheritParams F_habitats
#' @return a [list]
#' @export
F_habitats.forced <- function(t, pars){
  with(pars$QHABpar,{
    return(habitats_scale*habitats_t(t, pars))
  })
}

#' @title Set up forcing for the availability of habitats
#' @param pars a [list]
#' @param habitats_Opts a list, overwrites default values
#' @param habitats_scale a non-negative numeric value to scale the mean availability of habitats
#' @param habitats_t a function that specifies the temporal pattern for the availability of habitats
#' @return none
#' @export
setup_habitats_forced <- function(pars, habitats_Opts = list(), habitats_scale=0, habitats_t=NULL) {with(habitats_Opts,{

  QHABpar <- list()
  class(QHABpar) <- 'foi'
  pars$QHABpar <- QHABpar

  pars$QHABpar$habitats_scale = habitats_scale
  if(is.null(habitats_t)) habitats_t = function(t, pars){1}
  pars$QHABpar$habitats_t = habitats_t

  return(pars)
})}

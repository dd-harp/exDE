# specialized methods for when resource availability does not change dynamically

#' @title Modify resources and resource availability
#' @description Implements [Resources] for static models
#' @inheritParams Resources
#' @return none
#' @export
Resources.static <- function(t, y, pars) {
  return(pars)
}

#' @title Update blood feeding for static models
#' @description Implements [update_BloodFeeding] for a dynamic model
#' @inheritParams update_BloodFeeding
#' @return a [list]
#' @export
update_BloodFeeding.static <- function(t, y, pars){
  return(pars)
}

#' @title Update blood feeding for static models
#' @description Implements [update_EggLaying] for a dynamic model
#' @inheritParams update_EggLaying
#' @return a [list]
#' @export
update_EggLaying.static <- function(t, y, pars){
  return(pars)
}

#' @title Update host availability for static models
#' @description Implements [update_W] when pars$W doesn't change
#' @inheritParams update_W
#' @return a [list]
#' @export
update_W.static <- function(t, y, pars){
  return(pars)
}

#' @title Update other host availability for static models
#' @description Implements [update_O] when pars$O doesn't change
#' @inheritParams update_O
#' @return a [list]
#' @export
update_O.static <- function(t, pars){
  return(pars)
}

#' @title Update sugar availability for static models
#' @description Implements [update_S] when pars$S doesn't change
#' @inheritParams update_S
#' @return a [list]
#' @export
update_S.static <- function(t, pars){
  return(pars)
}

#' @title Set up parameters for a static model for resource availability
#' @param pars a [list]
#' @param RAopts a set of options that would override the defaults
#' @param other the availability of other blood hosts
#' @param sugar the availability of sugar
#' @return none
#' @export
setup_resources_static <- function(pars, RAopts, other=0, sugar=0){with(RAopts,{
  RApar <- list()
  class(RApar) <- 'static'
  pars$RApar <- RApar

  pars$OtherBloodHosts <-  other
  class(pars$OtherBloodHosts) <- 'static'

  pars$W <- compute_W(0, 0, pars)
  class(pars$W) <- 'static'

  pars$local_frac <- compute_local_frac(pars)

  pars$B <- compute_B(pars)
  class(pars$B) <- 'static'

  pars$Q <- compute_Q(pars)
  class(pars$Q) <- 'static'

  pars$S <- sugar
  class(pars$S) <- 'static'

  return(pars)
})}


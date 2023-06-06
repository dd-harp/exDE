# specialized methods for models with dynamically changing resource availability

#' @title Modify resources and resource availability dynamically
#' @description Implements [Resources] for the null model of resources (do nothing)
#' @inheritParams Resources
#' @return none
#' @export
Resources.dynamic <- function(t, y, pars) {
  pars = update_BloodFeeding(pars)
  pars = update_EggLaying(t, y, pars)
  pars = update_S(t, pars)
  return(pars)
}

#' @title Update total availability of blood hosts for dynamical models
#' @description Implements [update_BloodFeeding] for a dynamic model
#' @inheritParams update_BloodFeeding
#' @return a [list]
#' @export
update_BloodFeeding.dynamic <- function(t, y, pars){
  pars = update_W(t, y, pars)
  pars = update_O(t, pars)
  pars$B = compute_B(pars)
  pars$local_frac = compute_local_frac(pars)
  return(pars)
}

#' @title Update host availability dynamically
#' @description Implements [update_W] for a static model
#' @inheritParams update_W
#' @return a [list]
#' @export
update_W.dynamic<- function(t, y, pars){
  pars$W = compute_W(t, y, pars)
  return(pars)
}

#' @title Update other host availability dynamically
#' @description Implements [update_O] for a static model
#' @inheritParams update_O
#' @return a [list]
#' @export
update_O.dynamic<- function(t, pars){
  pars$OtherBloodHosts = F_otherblood(t, pars)
  return(pars)
}

#' @title Update the total availability of aquatic habitats and egg laying
#' @description Implements [update_EggLaying] for a dynamic model
#' @inheritParams update_EggLaying
#' @return a [list]
#' @export
update_EggLaying.dynamic <- function(t, y, pars){
  pars = F_habitats(t, pars)
  pars$Q = compute_Q(pars)
  pars$calU = make_calU(pars$calN, pars$searchQ)
  return(pars)
}

#' @title Dynamically update the total availability of sugar
#' @description Implements [update_S] for a dynamic model
#' @inheritParams update_S
#' @return a [list]
#' @export
update_S.dynamic <- function(t, pars){
  pars$S = F_sugar(t, pars)
  return(pars)
}

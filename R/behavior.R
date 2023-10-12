# human and mosquito behaviors

#' @title Methods for dynamic human and mosquito behaviors
#' @description This method dispatches on the type of `pars$BEHAVIOR`.
#' @param t current simulation time
#' @param y state variables
#' @param pars a [list]
#' @return [list]
#' @export
Behavior <- function(t, y, pars) {
  UseMethod("Behavior", pars$BEHAVIOR)
}

#' @title Methods for dynamic human and mosquito behaviors
#' @description Implements [Behavior] for the null model (no changes)
#' @inheritParams Behavior
#' @return [list]
#' @export
Behavior.null <- function(t, y, pars) {
 return(pars)
}

#' @title Make parameters for the null model for resource availability (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_behavior_null<- function(pars) {
  BEHAVIOR <- list()
  class(BEHAVIOR) <- 'null'
  pars$BEHAVIOR <- BEHAVIOR
  return(pars)
}

#' @title Setup behavior
#' @param pars a [list]
#' @return [list]
#' @export
setup_behavior <- function(pars) {
  UseMethod("setup_behavior", pars$BEHAVIOR)
}

#' @title Setup behavior
#' @param pars a [list]
#' @return [list]
#' @export
setup_behavior.null <- function(pars) {
  setup_behavior_forced(pars)
}

#' @title Setup behavior
#' @param pars a [list]
#' @return [list]
#' @export
setup_behavior.forced<- function(pars) {pars}

#' @title Methods for dynamic human and mosquito behaviors
#' @description Implements [Behavior] for the forced model (no changes)
#' @inheritParams Behavior
#' @return [list]
#' @export
Behavior.forced <- function(t, y, pars) {
  pars = UseBedNet(t, y, pars)
  pars = CareSeeking(t, y, pars)
#  pars = Protect(t, y, pars)
#  pars = Mobility(t, y, pars)
#  pars = MozySearch(t, y, pars)
 return(pars)
}

#' @title Make parameters for the forced model for resource availability (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_behavior_forced<- function(pars) {
  BEHAVIOR <- list()
  class(BEHAVIOR) <- 'forced'
  pars$BEHAVIOR <- BEHAVIOR
  pars = setup_care_seeking_null(pars)
#  pars = setup_protect_null(t, y, pars)
#  pars = setup_mobility_null(t, y, pars)
#  pars = setup_mozy_search_null(t, y, pars)
  return(pars)
}

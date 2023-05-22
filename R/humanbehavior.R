#' @title Simulate Human Behavior
#' @description This method dispatches on the type of `pars$HBpar`.
#' @param t current simulation time
#' @param y state variables
#' @param pars a [list]
#' @return none
#' @export
HumanBehavior <- function(t, y, pars) {
  UseMethod("HumanBehavior", pars$HBpar)
}

#' @title Simulate no human behaviors
#' @description Implements [HumanBehavior] for the null model (do nothing)
#' @inheritParams HumanBehavior
#' @return none
#' @export
HumanBehavior.null <- function(t, y, pars) {
  return(pars)
}

#' @title Make parameters for the null model for resource availability (do nothing)
#' @param pars a [list]
#' @return none
#' @export
setup_humanbehavior_null<- function(pars) {
  HBpar <- list()
  class(HBpar) <- 'null'
  pars$HBpar <- HBpar
  return(pars)
}

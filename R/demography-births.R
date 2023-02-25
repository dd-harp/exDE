
#' @title A function that computes the human population birth rate
#' @description This method dispatches on the type of `pars$Hpar$birthF`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return see help pages for specific methods
#' @export
F_births <- function(t, y, pars){
  UseMethod("Births", pars$Hpar$birthF)
}

#' @title A function that computes the human population birth rate
#' @description Implements [F_births] for the constant births model
#' @inheritParams F_births
#' @return a [numeric] vector of length `nStrata`
#' @export
F_births.0 = function(t, y, pars){with(pars$Hpar,{
  birthrate
})}

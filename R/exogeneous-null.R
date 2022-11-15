# specialized methods for the null model of exogenous forcing

#' @title Modify parameters due to exogenous forcing
#' @description Implements [ExogenousForcing] for the null model of exogeneous forcing (do nothing)
#' @inheritParams ExogenousForcing
#' @return nothing
#' @export
ExogenousForcing.null <- function(t, y, pars) {}

#' @title Make parameters for the null model of exogenous forcing (do nothing)
#' @param pars an [environment]
#' @return none
#' @export
make_parameters_exogenous_null <- function(pars) {
  stopifnot(is.environment(pars))
  EXpar <- list()
  class(EXpar) <- 'null'
  pars$EXpar <- EXpar
}

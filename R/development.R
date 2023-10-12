# generic methods for developments

#' @title Set up developments
#' @description This method dispatches on the type of `pars$DEVELOPMENT`.
#' @param t current simulation time
#' @param pars a [list]
#' @return [list]
#' @export
Development <- function(t, pars) {
  UseMethod("Development", pars$DEVELOPMENT)
}

#' @title Set up developments
#' @description Implements [Development] for the null model (do nothing)
#' @inheritParams Development
#' @return [list]
#' @export
Development.null <- function(t, pars) {
  return(pars)
}

#' @title Set up the null model for developments (do nothing)
#' @param pars a [list]
#' @return [list]
#' @export
setup_development_null <- function(pars) {
  DEVELOPMENT <- list()
  class(DEVELOPMENT) <- 'null'
  pars$DEVELOPMENT <- DEVELOPMENT
  return(pars)
}

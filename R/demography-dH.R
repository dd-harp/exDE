
#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] when `y` is static
#' @inheritParams dHdt
#' @return a [numeric] vector of 0s
#' @export
dHdt.zero <- function(t, y, pars, i){
  0*y
}

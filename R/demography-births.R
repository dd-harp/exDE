
#' @title Derivatives of demographic changes in human populations
#' @description Implements [Births] for static population models
#' @inheritParams Births
#' @return a [numeric] vector of zeros
#' @export
Births.zero <- function(t, y, pars){
  0*y
}

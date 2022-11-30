# dynamic model of demography

#' @title Size of human population denominators
#' @description Implements [F_H] for dynamic models.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.dynamic <- function(t, y, pars) {
  y[pars$H_ix]
}

#' @title Derivatives of demographic changes in human populations
#' @description Implements [dHdt] for dynamic models.
#' @inheritParams dHdt
#' @return a [numeric] vector of length `nStrata`
#' @export
dHdt.dynamic <- function(pars, ...) {
  n <- ...length()
  stopifnot(n %% 2 == 0) # even number of args
  DYs <- lapply(X = 1:(n/2), FUN = function(x) {
    D <- ...elt(x*2-1)
    stopifnot(is.matrix(D))
    Y <- ...elt(x*2)
    D %*% Y
  })
  rowSums(do.call(what = 'cbind', args = DYs))
}

#' @title Add indices for human population denominators to parameter list
#' @description Implements [make_index_H] for dynamic models.
#' @inheritParams make_index_H
#' @return none
#' @importFrom utils tail
#' @export
make_index_H.dynamic <- function(pars) {
  pars$H_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$H_ix, 1)
}

#' @title Make parameters for forced (trace) human demography model
#' @param pars an [environment]
#' @param b a list containing vectors of birth rates for each state in \eqn{\mathcal{X}}; each vector should be of length `nStrata`
#' @param d a list containing vectors of death rates for each state in \eqn{\mathcal{X}}; each vector should be of length `nStrata`
#' @param m a list containing vectors of ageing rates for each state in \eqn{\mathcal{X}}; each vector should be of length `nStrata-1`
#' @return none
#' @export
make_parameters_demography_dynamic <- function(pars, b, d, m) {
  stopifnot(is.environment(pars))
  stopifnot(is.list(b), is.list(d), is.list(m))
  stopifnot(all(vapply(b, function(x) {length(x) == pars$nStrata}, logical(1))))
  stopifnot(all(vapply(d, function(x) {length(x) == pars$nStrata}, logical(1))))
  stopifnot(all(vapply(m, function(x) {length(x) == pars$nStrata - 1}, logical(1))))
  Hpar <- list()
  class(Hpar) <- c('dynamic')
  Hpar$b <- b
  Hpar$d <- d
  Hpar$m <- m
  pars$Hpar <- Hpar
}

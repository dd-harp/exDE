
#' @title Make base parameters
#' @param solve_as, either "ode" or "dde"
#' @return a [list]
#' @export
make_parameters_xde = function(solve_as='ode'){
  pars = list()
  xde <- 'ode'
  class(xde) <- xde
  pars$xde = xde
  return(pars)
}

#' @title Set indices for generalized spatial model
#' @param pars an [environment]
#' @return none
#' @export
make_indices <- function(pars) {
  pars$max_ix <- 0
  if ('Lpar' %in% names(pars)) {
    pars = make_indices_L(pars)
  }
  if ('MYZpar' %in% names(pars)) {
    pars = make_indices_MYZ(pars)
  }
  if ('Xpar' %in% names(pars)) {
    pars = make_indices_X(pars)
  }
  if ('Hpar' %in% names(pars)) {
    pars = make_indices_H(pars)
  }
  return(pars)
}

#' @title Get the initial values as a vector
#' @param pars an [environment]
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_inits <- function(pars){
  if ('Lpar' %in% names(pars)) {
    Li = get_inits_L(pars)
  } else {Li = numeric(0)}
  if ('MYZpar' %in% names(pars)) {
    MYZi = get_inits_MYZ(pars)
  } else {MYZi = numeric(0)}
  if ('Xpar' %in% names(pars)) {
    Xi = get_inits_X(pars)
  } else {Xi = numeric(0)}
  if ('Hpar' %in% names(pars)) {
    Hi = get_inits_H(pars)
  } else {Hi = numeric(0)}
  y = c(Li, MYZi, Xi, Hi)
  #class(y) <- "dynamic"
  return(y)
}

#' @title Invert a diagonal matrix
#' @description Invert a diagonal matrix which is passed as a vector. If any
#' elements are zero, set them to one.
#' @param x a [numeric] vector
#' @return a diagonal [matrix]
#' @export
diag_inverse <- function(x) {
  x <- as.vector(x)
  ix <- which(x == 0)
  if (length(ix) > 0) {
    x[ix] <- 1
  }
  return(diag(x = 1/x, nrow = length(x), ncol = length(x), names = FALSE))
}

#' @title Check if two numeric values are approximately equal
#' @param a a [numeric] object
#' @param b a [numeric] object
#' @param tol the numeric tolerance
#' @return a [logical] value
#' @export
approx_equal <- function(a, b, tol = sqrt(.Machine$double.eps)) {
  abs(a - b) < tol
}


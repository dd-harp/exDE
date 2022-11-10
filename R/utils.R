#' @title Set indices for generalized spatial model
#' @param pars a [list]
#' @return modified parameters [list]
#' @export
make_indices <- function(pars) {
  pars$max_ix <- 0
  if ('Lpar' %in% names(pars)) {
    pars <- make_index_L(pars)
  }
  if ('MYZpar' %in% names(pars)) {
    pars <- make_index_MYZ(pars)
  }
  if ('Xpar' %in% names(pars)) {
    pars <- make_index_X(pars)
  }
  return(pars)
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
#' @return a logical value
#' @export
approx_equal <- function(a, b, tol = sqrt(.Machine$double.eps)) {
  abs(a - b) < tol
}

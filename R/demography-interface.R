# generic methods for demography (nested within human; \cal{H} in \cal{X})

#' @title Size of human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H <- function(t, y, pars) {
  UseMethod("F_H", pars$Hpar)
}

#' @title Size of lagged human population denominators
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param lag duration of lag `t-lag`
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H_lag <- function(t, y, pars, lag) {
  UseMethod("F_H_lag", pars$Hpar)
}

#' @title Derivatives of demographic changes in human populations
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @param ... additional arguments which must be of the form `D`, `Y`, etc. `D`
#' is a `nStrata` by `nStrata` demography matrix and `Y` is a state vector. There
#' should be at least one pair of these passed to the function.
#' @return see help pages for specific methods
#' @export
dHdt <- function(pars, ...) {
  UseMethod("dHdt", pars$Hpar)
}

#' @title Add indices for human population denominators to parameter list
#' @description This method dispatches on the type of `pars$Hpar`.
#' @param pars an [environment]
#' @return none
#' @export
make_index_H <- function(pars) {
  UseMethod("make_index_H", pars$Hpar)
}

#' @title Make the demography matrix \eqn{\mathcal{D}}
#' @description Make the demography matrix \eqn{\mathcal{D}}, which multiplies a
#' state vector `Y` on the left to return the derivatives of the state with respect
#' to demographic dynamics. Please note that while lengths of arguments are given
#' below, the function does no argument checking as it is calculated internally
#' assuming correct inputs.
#' @param d a vector of death rates (should have length equal to `nStrata`)
#' @param m a vector of ageing rates (should have length equal to `nStrata-1`)
#' @param b either `NULL` or a vector of birth rates (in which case should have length equal to `nStrata`)
#' @return a square [numeric] matrix of dimensions equal to the length of the input
#' argument `d`
#' @export
make_calD <- function(d, m, b = NULL) {
  calD <- matrix(data = 0, nrow = length(d), ncol = length(d))
  diag(calD) <- -(d + c(m, 0))
  # add birth rates if they exist
  if (!is.null(b)) {
    calD[1, ] <- calD[1, ] + b
  }
  # set subdiagonal for m
  calD[row(calD) == col(calD) + 1] <- m
  return(calD)
}

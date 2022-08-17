# specialized methods for the aquatic mosquito trace model

#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the trace (forced emergence) model.
#' @inheritParams F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.trace <- function(t, y, pars) {
  pars$Lpar$Lambda
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the trace (forced emergence) model.
#' @inheritParams dLdt
#' @return a [numeric] vector
#' @export
dLdt.trace <- function(t, y, pars, eta) {
  numeric(0)
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_index_L] for trace (forced emergence) model.
#' @inheritParams make_index_L
#' @return the modified parameter [list]
#' @export
make_index_L.trace <- function(pars) {
  pars$L_ix <- integer(0)
  return(pars)
}

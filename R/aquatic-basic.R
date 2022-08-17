# specialized methods for the aquatic mosquito basic competition model

#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the basic competition model.
#' @inheritParams F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.basic <- function(t, y, pars) {
  L <- y[pars$L_ix]
  pars$Lpar$psi*L
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the basic competition model.
#' @inheritParams dLdt
#' @return a [numeric] vector
#' @export
dLdt.basic <- function(t, y, pars, eta) {
  L <- y[pars$L_ix]
  with(pars$Lpar, {
    dL = eta - (psi + phi + (theta*L))*L
    return(dL)
  })
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_index_L] for basic competition model.
#' @inheritParams make_index_L
#' @return the modified parameter [list]
#' @importFrom utils tail
#' @export
make_index_L.basic <- function(pars) {
  pars$L_ix <- seq(from = pars$max_ix+1, length.out = pars$nHabitats)
  pars$max_ix <- tail(pars$L_ix, 1)
  return(pars)
}

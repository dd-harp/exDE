# specialized methods for the aquatic mosquito basic competition model

#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the basic competition model.
#' @inheritParams F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.basic <- function(t, y, pars) {
  with(pars$Lpar,{
    L <- y[L_ix]
    psi*L
  })
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the basic competition model.
#' @inheritParams dLdt
#' @return a [numeric] vector
#' @export
dLdt.basic <- function(t, y, pars, eta) {
  with(pars$Lpar, {
    L <- y[L_ix]
    dL = eta - (psi + phi + (theta*L))*L
    return(dL)
  })
}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_L] for basic competition model.
#' @inheritParams make_indices_L
#' @return none
#' @importFrom utils tail
#' @export
make_indices_L.basic <- function(pars) {
  pars$Lpar$L_ix <- seq(from = pars$max_ix+1, length.out = pars$nHabitats)
  pars$max_ix <- tail(pars$Lpar$L_ix, 1)
  return(pars)
}


#' @title Make parameters for basic competition aquatic mosquito model
#' @param pars an [list]
#' @param psi maturation rates for each aquatic habitat
#' @param phi density-independent mortality rates for each aquatic habitat
#' @param theta density-dependent mortality terms for each aquatic habitat
#' @return a [list] with Lpar added
#' @export
make_parameters_L_basic <- function(pars, psi, phi, theta) {
  stopifnot(is.numeric(psi), is.numeric(phi), is.numeric(theta))
  Lpar <- list()
  class(Lpar) <- 'basic'
  Lpar$psi <- psi
  Lpar$phi <- phi
  Lpar$theta <- theta
  pars$Lpar <- Lpar
  return(pars)
}

#' @title Make inits for basic competition aquatic mosquito model
#' @param pars an [list]
#' @param L0 initial conditions
#' @return a [list] with Linits added
#' @export
make_inits_L_basic <- function(pars, L0){
  stopifnot(is.numeric(L0))
  pars$Linits <- list(L0=L0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_L] for the GeRM model.
#' @inheritParams get_inits_L
#' @return none
#' @export
get_inits_L.basic <- function(pars){
  pars$Linits$L0
}


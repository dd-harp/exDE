# Compute the Transmission Terms

#' @title Compute the biting distribution matrix, beta
#' @description This method dispatches on the type of `pars$xde`
#' @param t current simulation time
#' @param y state vector
#' @param pars, a [list]
#' @return pars, a [list]
#' @export
F_beta <- function(t, y, pars){
  H <- F_H(t, y, pars)
  beta = compute_beta(H, pars$Hpar$wf, pars$Hpar$Psi)
  return(beta)
}

#' @title Entomological inoculation rate on human strata
#' @description This method dispatches on the type of `pars$Xpar$xde`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param beta, a [matrix]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars, beta) {
  fqZ <- with(pars$MYZpar,f*q)*F_Z(t, y, pars)
  fqZ <- fqZ_local(fqZ, pars)
  as.vector(beta %*% fqZ)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar$xde`.
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param beta, a [matrix]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa <- function(t, y, pars, beta) {
  x <- F_x(t, y, pars)
  kappa = t(beta) %*% x
  kappa = kappa_local(kappa, pars)
  as.vector(kappa)
}

#' @title Compute beta, the biting distribution matrix
#' @param H the human population size
#' @param wf the blood feeding search weights
#' @param Psi (time at risk), a [matrix]  dimensions `nPatches` by `nStrata`
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
compute_beta = function(H, wf, Psi){
  W <- as.vector(Psi %*% (wf*H))
  beta <- diag(wf, length(H)) %*% t(Psi) %*% diag(1/W, dim(Psi)[1])
  return(beta)
}

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
  beta = compute_beta(H, pars$Hpar$wts_f, pars$Hpar$TaR)
  return(beta)
}

#' @title Entomological inoculation rate on human strata
#' @description Compute the daily EIR for all the strata at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param beta, a [matrix]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars, beta) {
  fqZ <- F_fqZ(t, y, pars) * pars$local_frac
  as.vector(beta %*% fqZ)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Compute the net infectiousness of humans in each patch at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param beta, a [matrix]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa <- function(t, y, pars, beta) {
  kappa = as.vector(t(beta) %*% F_X(t, y, pars))
  with(pars, return(local_frac*kappa + (1-local_frac)*x_visitors))
}

#' @title Compute beta, the biting distribution matrix
#' @param H the human population size
#' @param wts_f the blood feeding search weights
#' @param TaR (time at risk), a [matrix]  dimensions `nPatches` by `nStrata`
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
compute_beta = function(H, wts_f, TaR){
  W <- as.vector(TaR %*% (wts_f*H))
  beta <- diag(wts_f, length(H)) %*% t(TaR) %*% diag(1/W, dim(TaR)[1])
  return(beta)
}

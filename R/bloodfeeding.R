# utilities to configure blood feeding

#' @title Entomological inoculation rate on human strata
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param MosyBehavior values returned by [exDE::MosquitoBehavior], potentially modified by control [exDE::VectorControl]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR <- function(t, y, pars, MosyBehavior) {
  Z <- F_Z(t, y, pars)
  f <- MosyBehavior$f[1]
  q <- MosyBehavior$q[1]
  beta <- F_beta(t, y, pars)
  as.vector(beta %*% diag(f*q, nrow = pars$nPatches) %*% Z)
}

#' @title Entomological inoculation rate on human strata at a lag
#' @param t lagged simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param MosyBehavior values returned by [exDE::MosquitoBehavior], potentially modified by control [exDE::VectorControl]
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR_lag <- function(t, y, pars, MosyBehavior){
  Z <- F_Z_lag(t, y, pars)
  f <- MosyBehavior$f[1]
  q <- MosyBehavior$q[1]
  beta <- F_beta_lag(t, y, pars)
  as.vector(beta %*% diag(f*q, nrow = pars$nPatches) %*% Z)
}

#' @title Biting distribution matrix
#' @description This method dispatches on the type of `pars$BFpar`
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta <- function(t, y, pars) {
  UseMethod("F_beta", pars$BFpar)
}

#' @title Lagged biting distribution matrix
#' @description This method dispatches on the type of `pars$BFpar`
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param lag is the lag value
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta_lag <- function(t, y, pars, lag) {
  UseMethod("F_beta_lag", pars$BFpar)
}

#' @title Biting distribution matrix
#' @description Implements [F_beta] when the availability of blood hosts is static
#' @inheritParams F_beta
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta.static <- function(t, y, pars) {
  return(pars$BFpar$beta)
}

#' @title Biting distribution matrix
#' @description Implements [F_beta_lag] when the availability of blood hosts is static
#' @inheritParams F_beta_lag
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta_lag.static <- function(t, y, pars, lag) {
  return(pars$BFpar$beta)
}

#' @title Biting distribution matrix
#' @description Implements [F_beta] for static denominators
#' @inheritParams F_beta
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta.dynamic <- function(t, y, pars) {
  H <- F_H(t, y, pars)
  W <- as.vector(pars$Hpar$Psi %*% (pars$Hpar$wf * H))
  return(
    diag(pars$Hpar$wf, pars$nStrata) %*% t(pars$Hpar$Psi) %*% diag(1/W, pars$nPatches)
  )
}

#' @title Lagged biting distribution matrix
#' @description Implements [F_beta_lag] for dynamic denominators
#' @inheritParams F_beta_lag
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
F_beta_lag.dynamic <- function(t, y, pars, lag) {
  H <- F_H_lag(t, y, pars, lag)
  W <- as.vector(pars$Hpar$Psi %*% (pars$Hpar$wf * H))
  return(
    diag(pars$Hpar$wf, pars$nStrata) %*% t(pars$Hpar$Psi) %*% diag(1/W, pars$nPatches)
  )
}


#' @title Sets up BFpars for a model with static availability
#' @param params is a list
#' @return a [list] with BFpar added
#' @export
make_parameters_BF_static <- function(params){with(params$Hpar,{
  BFpar <- list()
  class(BFpar) <- c('static')
  H <- F_H(0,0,params)
  TaR <- TimeSpent
  W <- as.vector(TaR %*% (searchWtsH*H))
  BFpar$beta = diag(searchWtsH) %*% t(TaR) %*% diag(1/W)
  params$BFpar = BFpar
  return(params)
})}

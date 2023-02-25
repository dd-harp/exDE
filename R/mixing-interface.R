# generic methods to compute mixing matrices, beta and beta_lag

#' @title Compute beta, the biting distribution matrix
#' @description This method dispatches on the type of `pars$Hpar`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return pars, a [list]
#' @export
make_beta <- function(t, y, pars) {
  UseMethod("make_beta", pars$Hpar)
}

#' @title Compute beta_lag, the lagged biting distribution matrix
#' @description This method dispatches on the type of `pars$Hpar`
#' @param t current simulation time
#' @param y state vector
#' @param pars an [environment]
#' @param lag is the lag value
#' @return a [matrix] of dimensions `nStrata` by `nPatches`
#' @export
make_beta_lag <- function(t, y, pars, lag) {
  UseMethod("make_beta_lag", pars$Hpar)
}

#' @title Compute beta, the biting distribution matrix
#' @description This method dispatches on the type of `pars$Hpar`
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

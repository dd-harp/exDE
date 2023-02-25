# specialized methods to compute beta with dynamic denominators

#' @title Biting distribution matrix
#' @description Implements [make_beta] for dynamic denominators
#' @description Compute beta, a [matrix] of dimensions `nStrata` by `nPatches`
#' @inheritParams make_beta
#' @return pars, a [list]
#' @export
make_beta.dynamic <- function(t, y, pars) {
  H <- F_H(t, y, pars)
  pars$beta = compute_beta(H, pars$Hpar$wf, pars$Hpar$Psi)
  return(pars)
}

#' @title Lagged biting distribution matrix
#' @description Implements [make_beta] for dynamic denominators
#' @description Compute beta, a [matrix] of dimensions `nStrata` by `nPatches`
#' @inheritParams make_beta_lag
#' @return pars, a [list]
#' @export
make_beta_lag.dynamic <- function(t, y, pars, lag) {
  H_lag <- F_H_lag(t, y, pars, lag)
  pars$beta_lag = compute_beta(H_lag, pars$Hpar$wf, pars$Hpar$Psi)
  return(pars)
}

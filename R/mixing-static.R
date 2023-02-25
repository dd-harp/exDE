# specialized methods to compute beta with static denominators

#' @title Biting distribution matrix
#' @description Implements [make_beta] for static denominators
#' @inheritParams make_beta
#' @return pars, a [list]
#' @export
make_beta.static <- function(t, y, pars) {
  return(pars)
}

#' @title Biting distribution matrix
#' @description Implements [make_beta_lag] for static denominators
#' @inheritParams make_beta_lag
#' @return pars, a [list]
#' @export
make_beta_lag.static <- function(t, y, pars, lag){
  return(pars)
}

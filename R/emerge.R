#' @title Compute emerging adults
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
#' @export
Emergence = function(t, y, pars){
  for(s in 1:pars$nVectors)
    pars$Lambda[[s]] = pars$calN %*% F_alpha(t, y, pars,s)
  return(pars)
}

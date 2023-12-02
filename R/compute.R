#' @title Compute the EIR
#' @description Using the output of deSolve,
#' compute the EIR at every point in time
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_EIR <- function(deout, pars) {
  ix = 1:length(deout[,1])
  eir = sapply(ix, compute_EIR_i, deout=deout, pars=pars)
  return(eir)
}

#' @title Compute the EIR for the ith point in time
#' @description Using the output of deSolve,
#' compute the EIR at one point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric]
#' @export
compute_EIR_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  compute_EIR_ty(t, y, pars)
}

#' @title Compute the EIR
#' @description Using the output of deSolve,
#' compute the EIR at one point in time
#' @param t [numeric]
#' @param y [vector] the state variables
#' @param pars a [list]
#' @return [numeric]
#' @export
compute_EIR_ty <- function(t, y, pars) {
  beta <- F_beta(t, y, pars)
  eir = F_EIR(t, y, pars, beta)
  return(eir)
}

#' @title Compute the fqZ
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the fqZ at every point in time
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [vector]
#' @export
compute_fqZ <- function(deout, pars) {
  ix = 1:length(deout[,1])
  fqZ = sapply(ix, compute_fqZ_i, deout=deout, pars=pars)
  return(fqZ)
}

#' @title Compute the fqZ for the ith
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the fqZ at one point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric]
#' @export
compute_fqZ_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  fqZ = F_fqZ(t, y, pars)
  return(fqZ)
}

#' @title Compute the NI
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the NI for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric] containing the NI
#' @export
compute_NI <- function(deout, pars) {
  ix = 1:length(deout[,1])
  NI = sapply(ix, compute_NI_i, deout=deout, pars=pars)
  return(NI)
}

#' @title Compute NI once
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the NI at a point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return a [list] containing the NI
#' @export
compute_NI_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  X <- F_X(t, y, pars)
  H = F_H(t, y, pars)
  NI = X/H
  return(NI)
}

#' @title Compute kappa
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return [numeric] containing the kappa
#' @export
compute_kappa <- function(deout, pars) {
  ix = 1:length(deout[,1])
  kappa = sapply(ix, compute_kappa_i, deout=deout, pars=pars)
  return(kappa)
}

#' @title Compute kappa once
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa at a point in time
#' @param i an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @return a [list] containing the kappa
#' @export
compute_kappa_i <- function(i, deout, pars) {
  t = deout[i,1]
  y = deout[i,-1]
  compute_kappa_ty(t,y,pars)
}

#' @title Compute kappa, ty
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa at a point in time
#' @param t [numeric]
#' @param y [vector] the state variables
#' @param pars a [list]
#' @return a [list] containing the kappa
#' @export
compute_kappa_ty <- function(t, y, pars) {
  beta <- F_beta(t, y, pars)
  kappa = F_kappa(t, y, pars, beta)
  return(kappa)
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve compute the dynamical terms every point in time. This method dispatches on the type of `pars$compute`.
#' @param varslist a [list] with variable names attached
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_terms <- function(varslist, deout, pars) {
  UseMethod("compute_terms", pars$compute)
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve.ode` or `xde_solve.dde`
#' @param varslist a [list] with variable names attached
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_terms.xde <- function(varslist, deout, pars) {
  eir = compute_EIR(deout, pars)
  ni = compute_NI(deout, pars)
  kappa = compute_kappa(deout, pars)
  fqZ = compute_fqZ(deout, pars)
  pr = F_pr(varslist, pars)
  return(list(eir=eir,pr=pr,ni=ni,kappa=kappa,fqZ=fqZ))
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve.cohort`
#' @param varslist a [list] with variable names attached
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_terms.cohort <- function(varslist, deout, pars) {
  eir = pars$F_eir(deout[,1], pars)
  ni = compute_NI(deout, pars)
  pr = F_pr(varslist, pars)
  return(list(eir=eir,pr=pr,ni=ni))
}


#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve.human`
#' @param varslist a [list] with variable names attached
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_terms.human<- function(varslist, deout, pars) {
  eir = compute_EIR(deout, pars)
  ni = compute_NI(deout, pars)
  fqZ = compute_fqZ(deout, pars)
  pr = F_pr(varslist, pars)
  return(list(eir=eir,pr=pr,ni=ni,fqZ=fqZ))
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve` that
#' don't have any relevant dynamical terms
#' @param varslist a [list] with variable names attached
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_terms.na <- function(varslist,deout, pars) {
 return(list())
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `root_solve`
#' @param varslist a [list] with variable names attached
#' @param y_eq a [matrix], the output of `rootSolve`
#' @param pars a [list]
#' @return [matrix]
#' @export
compute_terms_steady<- function(varslist, y_eq, pars) {
  eir <- compute_EIR_ty(0, y_eq, pars)
  kappa <- compute_kappa_ty(0, y_eq, pars)
  fqZ <- F_fqZ(0, y_eq, pars)
  ni <- F_X(0, y_eq, pars)/varslist$XH$H
  pr <- F_pr(varslist, pars)/varslist$XH$H
  return(list(eir=eir,pr=pr,kappa=kappa,fqZ=fqZ,ni=ni))
}

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

#' @title Compute the EIR
#' @description Using the output of deSolve,
#' compute the EIR at the i^th point in time
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
#' @description Using the time and a vector of variables,
#' compute the EIR
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

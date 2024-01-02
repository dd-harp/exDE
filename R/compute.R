#' @title Compute dynamical terms
#' @description Using the output of deSolve compute the dynamical terms every point in time. This method dispatches on the type of `pars$compute`.
#' @param varslist a [list] with variable names attached
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @param s the vector species index
#' @param i the host species index
#' @return [matrix]
#' @export
compute_terms <- function(varslist, deout, pars, s, i) {
  UseMethod("compute_terms", pars$compute)
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve.ode` or `xde_solve.dde`
#' @inheritParams compute_terms
#' @return [list]
#' @export
compute_terms.xde <- function(varslist, deout, pars, s, i) {
  time = deout[,1]
  eir = compute_EIR(deout, pars, s, i)
  ni = compute_NI(deout, pars, i)
  kappa = compute_kappa(deout, pars, i)
  fqZ = compute_fqZ(deout, pars, s)
  pr = F_pr(varslist, pars, i)
  return(list(time=time,eir=eir,pr=pr,ni=ni,kappa=kappa,fqZ=fqZ))
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve.cohort`
#' @inheritParams compute_terms
#' @return [list]
#' @export
compute_terms.cohort <- function(varslist, deout, pars, s, i) {
  time = deout[,1]
  d1 = length(time)
  eir = matrix(pars$F_eir(time, pars), d1, pars$nStrata)
  ni = compute_NI(deout, pars, i)
  pr = F_pr(varslist, pars, i)
  return(list(time=time,eir=eir,pr=pr,ni=ni))
}


#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve.human`
#' @inheritParams compute_terms
#' @return [list]
#' @export
compute_terms.human<- function(varslist, deout, pars, s, i) {
  time = deout[,1]
  eir = compute_EIR(deout, pars, s, i)
  ni = compute_NI(deout, pars)
  fqZ = compute_fqZ(deout, pars)
  pr = F_pr(varslist, pars)
  return(list(time=time,eir=eir,pr=pr,ni=ni,fqZ=fqZ))
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `xde_solve` that
#' don't have any relevant dynamical terms
#' @inheritParams compute_terms
#' @return [matrix]
#' @export
compute_terms.na <- function(varslist,deout, pars, s, i) {
 return(list())
}

#' @title Compute dynamical terms
#' @description Using the output of deSolve
#' compute the dynamical terms for the output of `root_solve`
#' @param varslist a [list] with variable names attached
#' @param y_eq a [matrix], the output of `rootSolve`
#' @param pars a [list]
#' @param s the vector species index
#' @param i the host species index
#' @return [matrix]
#' @export
compute_terms_steady<- function(varslist, y_eq, pars, s, i) {
  eir <- compute_EIR_ty(0, y_eq, pars, s, i)
  kappa <- compute_kappa_ty(0, y_eq, pars, i)
  fqZ <- F_fqZ(0, y_eq, pars, s)
  ni <- F_X(0, y_eq, pars, i)/varslist$XH$H
  pr <- F_pr(varslist, pars, i)
  return(list(eir=eir,pr=pr,kappa=kappa,fqZ=fqZ,ni=ni))
}

#' @title Compute the EIR
#' @description Using the output of deSolve,
#' compute the EIR at every point in time
#' @param deout a [matrix], the output of deSolve
#' @param pars a [list]
#' @param s the vector species index
#' @param i the host species index
#' @return [matrix]
#' @export
compute_EIR <- function(deout, pars, s, i) {
  d1 = length(deout[,1])
  eir = sapply(1:d1, compute_EIR_ix, deout=deout, pars=pars, s=s, i=i)
  eir = shapeIt(eir, d1, pars$nStrata)
  return(eir)
}

#' @title Compute the EIR
#' @description Using the output of deSolve,
#' compute the EIR at the i^th point in time
#' @param ix an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @param s the vector species index
#' @param i the host species index
#' @return [numeric]
#' @export
compute_EIR_ix <- function(ix, deout, pars, s, i) {
  t = deout[ix,1]
  y = deout[ix,-1]
  compute_EIR_ty(t, y, pars, s, i)
}

#' @title Compute the EIR
#' @description Using the time and a vector of variables,
#' compute the EIR
#' @param t [numeric]
#' @param y [vector] the state variables
#' @param pars a [list]
#' @param s the vector species index
#' @param i the host species index
#' @return [numeric]
#' @export
compute_EIR_ty <- function(t, y, pars, s, i) {
  beta <- F_beta(t, y, pars, i)
  eir = F_EIR(t, y, pars, beta, s)
  return(eir)
}

#' @title Compute the fqZ
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the fqZ at every point in time
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @param s the vector species index
#' @return [vector]
#' @export
compute_fqZ <- function(deout, pars, s) {
  d1 = length(deout[,1])
  fqZ = sapply(1:d1, compute_fqZ_ix, deout=deout, pars=pars, s=s)
  fqZ = shapeIt(fqZ, d1, pars$nPatches)
  return(fqZ)
}

#' @title Compute the fqZ for the ith
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the fqZ at one point in time
#' @param ix an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @param s the vector species index
#' @return [numeric]
#' @export
compute_fqZ_ix <- function(ix, deout, pars, s) {
  t = deout[ix,1]
  y = deout[ix,-1]
  fqZ = F_fqZ(t, y, pars, s)
  return(fqZ)
}

#' @title Compute the NI
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the NI for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @param i the host species index
#' @return [numeric] containing the NI
#' @export
compute_NI <- function(deout, pars, i) {
  d1 = length(deout[,1])
  NI = sapply(1:d1, compute_NI_ix, deout=deout, pars=pars, i=i)
  NI = shapeIt(NI, d1, pars$nStrata)
  return(NI)
}

#' @title Compute NI once
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the NI at a point in time
#' @param ix an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @param i the host species index
#' @return a [list] containing the NI
#' @export
compute_NI_ix <- function(ix, deout, pars, i) {
  t = deout[ix,1]
  y = deout[ix,-1]
  X <- F_X(t, y, pars, i)
  H = F_H(t, y, pars, i)
  NI = X/H
  return(NI)
}

#' @title Compute kappa
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa for each stratum
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @param i the host species index
#' @return [numeric] containing the kappa
#' @export
compute_kappa <- function(deout, pars, i) {
  d1 = length(deout[,1])
  kappa = sapply(1:d1, compute_kappa_ix, deout=deout, pars=pars, i=i)
  kappa = shapeIt(kappa, d1, pars$nPatches)
  return(kappa)
}

#' @title Compute kappa once
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa at a point in time
#' @param ix an [integer]
#' @param deout a matrix, the output of deSolve
#' @param pars a [list]
#' @param i the host species index
#' @return a [list] containing the kappa
#' @export
compute_kappa_ix <- function(ix, deout, pars, i) {
  t = deout[ix,1]
  y = deout[ix,-1]
  compute_kappa_ty(t,y,pars,i)
}

#' @title Compute kappa, ty
#' @description Using the output of [deSolve::ode] or [deSolve::dede],
#' compute the kappa at a point in time
#' @param t [numeric]
#' @param y [vector] the state variables
#' @param pars a [list]
#' @param i the host species index
#' @return a [list] containing the kappa
#' @export
compute_kappa_ty <- function(t, y, pars, i) {
  beta <- F_beta(t, y, pars,  i)
  kappa = F_kappa(t, y, pars, beta, i)
  return(kappa)
}

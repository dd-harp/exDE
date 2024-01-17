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

  eir = list()
  if(pars$nHosts>0) for(i in 1:pars$nHosts)
     eir[[i]] = matrix(0, nrow = length(time), ncol = pars$Hpar[[i]]$nStrata)

  kappa = list()
  fqZ = list()
  fqM = list()
  for(s in 1:pars$nVectors){
     kappa[[s]] = matrix(0, nrow = length(time), ncol = pars$nPatches)
     fqZ[[s]] = matrix(0, nrow = length(time), ncol = pars$nPatches)
     fqM[[s]] = matrix(0, nrow = length(time), ncol = pars$nPatches)
  }

  for (ix in 1:length(time)){
    yt = deout[ix,-1]
    pars = compute_vars_full(time[ix], yt, pars)

    for(i in 1:pars$nHosts)
      eir[[i]][ix,] = pars$EIR[[i]]


    for(s in 1:pars$nVectors){
      kappa[[s]][ix,] = pars$kappa[[i]]
      fqZ[[s]][ix,] = F_fqZ(time[ix], yt, pars, s)
      fqM[[s]][ix,] = F_fqM(time[ix], yt, pars, s)
    }
  }
  pr = list()
  ni = list()
  for(i in 1:pars$nHosts){
    ni[[i]] = compute_NI(deout, pars, i)
    pr[[i]] = F_pr(varslist, pars, i)
  }

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

  eir = matrix(pars$F_eir(time, pars), d1, pars$Hpar[[i]]$nStrata)
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

  eir = list()
  if(pars$nHosts>0) for(i in 1:pars$nHosts)
    eir[[i]] = matrix(0, nrow = length(time), ncol = pars$Hpar[[i]]$nStrata)

  fqZ = list()
  for(s in 1:pars$nVectors)
    fqZ[[s]] = matrix(0, nrow = length(time), ncol = pars$nPatches)

  for (ix in 1:length(time)){
    yt = deout[ix,-1]
    pars = compute_vars_full(time[ix], yt, pars)

    for(i in 1:pars$nHosts)
      eir[[i]][ix,] = pars$EIR[[i]]

    for(s in 1:pars$nVectors){
      kappa[[s]][ix,] = pars$kappa[[i]]
      fqZ[[s]][ix,] = F_fqZ(time[ix], yt, pars, s)
    }
  }

  pr = list()
  ni = list()
  for(i in 1:pars$nHosts){
    ni[[i]] = compute_NI(deout, pars, i)
    pr[[i]] = F_pr(varslist, pars, i)
  }

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
  pars = Transmission(0, y_eq, pars)
  eir = pars$EIR[[1]]
  fqZ <- F_fqZ(0, y_eq, pars, s)
  ni <- F_X(0, y_eq, pars, i)/varslist$XH[[i]]$H
  pr <- F_pr(varslist, pars, i)
  return(list(eir=eir,pr=pr,kappa=kappa,fqZ=fqZ,ni=ni))
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
  NI = shapeIt(NI, d1, pars$Hpar[[i]]$nStrata)
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


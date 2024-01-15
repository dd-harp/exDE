# a hybrid model tracking mean MoI for all and apparent infections

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the hybrid MoI model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @importFrom stats pexp
#' @export
F_X.hMoI <- function(t, y, pars, i) {
  with(pars$ix$X[[i]],{
  m1 = y[m1_ix]
  m2 = y[m2_ix]
  H <- F_H(t, y, pars, i)
  with(pars$Xpar[[i]],{
    x1 <- pexp(q = m1)
    x2 <- pexp(q = m2)
    X <- ((c2 * x2) + (c1 * (x1 - x2)))*H
    return(X)
  })
})}

#' @title Size of the human population
#' @description Implements [F_H] for the hybrid MoI model.
#' @inheritParams F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.hMoI <- function(t, y, pars, i) {
  pars$Hpar[[i]]$H
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the hMoI model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.hMoI<- function(varslist, pars, i) {
  pr = with(varslist$XH[[i]], 1-exp(-m1))
  return(pr)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the hMoI model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.hMoI <- function(y, pars,i) {
  with(pars$Xpar[[i]], b)
}


#' @title Derivatives for human population
#' @description Implements [dXdt] for the hybrid MoI model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.hMoI <- function(t, y, pars, i) {

  foi = pars$FoI[[i]]

  with(pars$ix$X[[i]],{
    m1 = y[m1_ix]
    m2 = y[m2_ix]

    with(pars$Xpar[[i]], {
      dm1dt <- foi - r1*m1
      dm2dt <- foi - r2*m2
      return(c(dm1dt, dm2dt))
    })
  })
}

#' @title Setup Xpar.hMoI
#' @description Implements [setup_Xpar] for the hMoI model
#' @inheritParams setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.hMoI = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_hMoI(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

#' @title Setup Xinits.hMoI
#' @description Implements [setup_Xinits] for the hMoI model
#' @inheritParams setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.hMoI = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = make_Xinits_hMoI(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}


#' @title Make parameters for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param nStrata is the number of human population strata
#' @param Xopts a [list] that overwrites default values
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c1 transmission probability (efficiency) from inapparent human infections to mosquito
#' @param c2 transmission probability (efficiency) from patent human infections to mosquito
#' @param r1 recovery rate from inapparent infections
#' @param r2 recovery rate from patent infections
#' @return none
#' @export
make_Xpar_hMoI = function(nStrata, Xopts=list(),
                          b=0.55, r1=1/180, r2 = 1/70,
                          c1=0.015, c2=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- "hMoI"

    Xpar$b = checkIt(b, nStrata)
    Xpar$c1 = checkIt(c1, nStrata)
    Xpar$c2 = checkIt(c2, nStrata)
    Xpar$r1 = checkIt(r1, nStrata)
    Xpar$r2 = checkIt(r2, nStrata)

    return(Xpar)
})}

#' @title Make inits for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param nStrata the number of population strata
#' @param Xopts a [list] that overwrites default values
#' @param m10 mean MoI among inapparent human infections
#' @param m20 mean MoI among patent human infections
#' @return none
#' @export
make_Xinits_hMoI = function(nStrata, Xopts = list(), m10=2, m20=1){with(Xopts,{
  m1 = checkIt(m10, nStrata)
  m2 = checkIt(m20, nStrata)
  return(list(m1=m1, m2=m2))
})}


#' @title Compute the HTC for the hMoI model
#' @description Implements [HTC] for the hMoI model with demography.
#' @inheritParams HTC
#' @return a [numeric] vector
#' @export
HTC.hMoI <- function(pars, i) {
  with(pars$Xpar[[i]],
       return(c2/r2 + c1*(1/r1 - 1/r2))
  )
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the hybrid MoI model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.hMoI <- function(pars, i) {with(pars,{
  m1_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(m1_ix, 1)

  m2_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(m2_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(m1_ix=m1_ix, m2_ix=m2_ix)
  return(pars)
})}

#' @title Make parameters for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param pars a [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c1 transmission probability (efficiency) from inapparent human infections to mosquito
#' @param c2 transmission probability (efficiency) from patent human infections to mosquito
#' @param r1 recovery rate from inapparent infections
#' @param r2 recovery rate from patent infections
#' @return none
#' @export
make_parameters_X_hMoI <- function(pars, b, c1, c2, r1, r2) {
  stopifnot(is.numeric(b), is.numeric(c1), is.numeric(c2), is.numeric(r1), is.numeric(r2))
  Xpar <- list()
  class(Xpar) <- c('hMoI')
  Xpar$b <- b
  Xpar$c1 <- c1
  Xpar$c2 <- c2
  Xpar$r1 <- r1
  Xpar$r2 <- r2
  pars$Xpar[[1]] <- Xpar
  return(pars)
}

#' @title Make inits for hybrid MoI human model
#' @description MoI stands for Multiplicity of Infection, and refers to malarial superinfection.
#' @param pars a [list]
#' @param m10 mean MoI among inapparent human infections
#' @param m20 mean MoI among patent human infections
#' @return none
#' @export
make_inits_X_hMoI <- function(pars, m10, m20) {
  stopifnot(is.numeric(m10), is.numeric(m20))
  pars$Xinits[[1]] = list(m1 = m10, m2 = m20)
  return(pars)
}

#' @title Update inits for hybrid MoI human model from a vector of states
#' @inheritParams update_inits_X
#' @return none
#' @export
update_inits_X.hMoI <- function(pars, y0,i) {
  with(pars$ix$X[[i]],{
    m1 = y0[m1_ix]
    m2 = y0[m2_ix]
    pars$Xinits[[i]] = make_Xinits_hMoI(pars$Hpar[[1]]$nStrata, m10=m1, m20=m2)
    return(pars)
})}

#' @title Parse the output of deSolve and return variables for the hMoI model
#' @description Implements [parse_deout_X] for the hMoI model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.hMoI <- function(deout, pars, i){
  H = pars$Hpar[[i]]$H
  with(pars$ix$X[[i]], {
    time = deout[,1]
    m1 = deout[,m1_ix+1]
    m2 = deout[,m2_ix+1]
  return(list(time=time, H=H,m1=m1,m2=m2))
})}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams get_inits_X
#' @return none
#' @export
get_inits_X.hMoI <- function(pars, i){with(pars$Xinits[[i]],{
  c(m1, m2)
})}

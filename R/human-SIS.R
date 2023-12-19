# specialized methods for the human SIS model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIS <- function(t, y, pars) {
  I = y[pars$ix$X$I_ix]
  X = with(pars$Xpar, c*I)
  return(X)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIS model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIS <- function(varslist, pars) {
  pr = with(varslist$XH, I/H)
  return(pr)
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIS model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIS <- function(y, pars) {
  with(pars$Xpar, b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model, no demography.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SISdX <- function(t, y, pars, FoI) {

  I <- y[pars$ix$X$I_ix]
  H <- F_H(t, y, pars)

  with(pars$Xpar, {
    dI <- FoI*(H - I) - r*I
    return(c(dI))
  })
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model with demography.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SISdXdH <- function(t, y, pars, FoI) {

  I <- y[pars$ix$X$I_ix]
  H <- F_H(t, y, pars)

  with(pars$Xpar, {

    dI <- FoI*(H - I) - r*I + dHdt(t, I, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    return(c(dI, dH))
  })
}

#' @title Setup Xpar.SIS
#' @description Implements [setup_X] for the SIS model
#' @inheritParams setup_X
#' @return a [list] vector
#' @export
setup_X.SIS = function(pars, Xname, Xopts=list()){

  pars$Xname = "SIS"
  pars = make_Xpar_SIS(pars, Xopts)
  pars = make_Xinits_SIS(pars, Xopts)

  return(pars)
}

#' @title Make parameters for SIS human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list]
#' @export
make_Xpar_SIS = function(pars, Xopts=list(),
                         b=0.55, r=1/180, c=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIS", "SISdX")

    Xpar$b = checkIt(b, pars$nStrata)
    Xpar$c = checkIt(c, pars$nStrata)
    Xpar$r = checkIt(r, pars$nStrata)

    pars$Xpar = Xpar
    return(pars)
  })}

#' @title Make initial values for the SIS human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] to overwrite defaults
#' @param I0 the initial values of the parameter X
#' @return a [list]
#' @export
make_Xinits_SIS = function(pars, Xopts = list(), I0=1){with(Xopts,{
  inits = list()
  inits$I0 = checkIt(I0, pars$nStrata)
  pars$Xinits = inits
  return(pars)
})}

#' @title Parse the output of deSolve and return variables for the SIS model
#' @description Implements [parse_deout_X] for the SIS model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.SIS <- function(deout, pars) {
  time = deout[,1]
  Hlist <- parse_deout_H(deout, pars)
  I = deout[,pars$ix$X$I_ix+1]
  with(Hlist,{
    return(list(time=time, I=I, H=H))
})}

#' @title Compute the HTC for the SIS model
#' @description Implements [HTC] for the SIS model with demography.
#' @inheritParams HTC
#' @return a [numeric] vector
#' @export
HTC.SIS <- function(pars) {
  with(pars$Xpar,
    return(c/r)
  )
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIS model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIS <- function(pars) {
  pars$ix$X$I_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$I_ix, 1)
  return(pars)
}

#' @title Make parameters for SIS human model
#' @param pars a [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list]
#' @export
make_parameters_X_SIS <- function(pars, b, c, r) {
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r))
  Xpar <- list()
  class(Xpar) <- c('SIS', 'SISdX')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  pars$Xpar <- Xpar
  return(pars)
}

#' @title Make inits for SIS human model
#' @param pars a [list]
#' @param I0 size of infected population in each strata
#' @return none
#' @export
make_inits_X_SIS <- function(pars, I0) {
  stopifnot(is.numeric(I0))
  pars$Xinits <- list(I0=I0)
  return(pars)
}

#' @title Update inits for the SIS human model from a vector of states
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X.SIS <- function(pars, y0) {
  I0 = y0[pars$ix$X$I_ix]
  make_Xinits_SIS(pars, list(), I0)
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X.SIS <- function(pars){
  pars$Xinits$I0
}

#' Plot the density of infected individuals for the SIS model
#'
#' @inheritParams xde_plot_X
#' @export
xde_plot_X.SIS = function(pars, clrs="black", llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH,
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X(vars$XH, pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SIS model
#'
#' @inheritParams xde_lines_X
#'
#' @export
xde_lines_X.SIS = function(XH, pars, clrs="black", llty=1){
  with(XH,{
    if(pars$nStrata==1) lines(time, X, col=clrs[1], lty = llty[1])
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, X[,i], col=clrs[i], lty = llty[i])
      }
    }
  })}

# specialized methods for the human SIS model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIS <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  X = with(pars$Xpar[[i]], c*I)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIS <- function(t, y, pars, i){
  S = y[pars$ix$X[[i]]$S_ix]
  I = y[pars$ix$X[[i]]$I_ix]
  return(S+I)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIS model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIS <- function(varslist, pars, i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIS model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIS <- function(y, pars, i) {
  with(pars$Xpar[[i]], b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model, no demography.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIS <- function(t, y, pars, FoI, i) {

  S <- y[pars$ix$X[[i]]$S_ix]
  I <- y[pars$ix$X[[i]]$I_ix]
  H <- F_H(t, y, pars, i)
  with(pars$Xpar[[i]], {
    dS <- Births(t, H, pars, i) - FoI*S + r*I + dHdt(t, S, pars, i)
    dI <- FoI*S - r*I + dHdt(t, I, pars, i)
    return(c(dS, dI))
  })
}

#' @title Setup Xpar.SIS
#' @description Implements [setup_Xpar] for the SIS model
#' @inheritParams setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SIS = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SIS(pars$nStrata, Xopts)
  return(pars)
}

#' @title Setup Xinits.SIS
#' @description Implements [setup_Xinits] for the SIS model
#' @inheritParams setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIS = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars,make_Xinits_SIS(nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}


#' @title Make parameters for SIS human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list]
#' @export
make_Xpar_SIS = function(nStrata, Xopts=list(),
                         b=0.55, r=1/180, c=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- "SIS"

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)

    return(Xpar)
  })}

#' @title Make initial values for the SIS human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial human population density
#' @param S0 the initial values of the parameter S
#' @param I0 the initial values of the parameter I
#' @return a [list]
#' @export
make_Xinits_SIS = function(nStrata, Xopts = list(), H0=NULL, S0=NULL, I0=1){with(Xopts,{
  if(is.null(S0)) S0 = H0 - I0
  stopifnot(is.numeric(S0))
  S = checkIt(S0, nStrata)
  I = checkIt(I0, nStrata)
  return(list(S=S, I=I))
})}

#' @title Parse the output of deSolve and return variables for the SIS model
#' @description Implements [parse_deout_X] for the SIS model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.SIS <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    I = deout[,I_ix+1]
    H = S+I
    return(list(time=time, S=S, I=I, H=H))
})}

#' @title Compute the HTC for the SIS model
#' @description Implements [HTC] for the SIS model with demography.
#' @inheritParams HTC
#' @return a [numeric] vector
#' @export
HTC.SIS <- function(pars, i) {
  with(pars$Xpar[[i]],
    return(c/r)
  )
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIS model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIS <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(I_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix)
  return(pars)
})}

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
  pars$Xpar[[1]] <- Xpar
  return(pars)
}

#' @title Make inits for SIS human model
#' @param pars a [list]
#' @param S0 size of infected population in each strata
#' @param I0 size of infected population in each strata
#' @return none
#' @export
make_inits_X_SIS <- function(pars, S0, I0) {
  stopifnot(is.numeric(S0))
  stopifnot(is.numeric(I0))
  pars$Xinits[[1]] <- list(S=S0, I=I0)
  return(pars)
}

#' @title Update inits for the SIS human model from a vector of states
#' @inheritParams update_inits_X
#' @return none
#' @export
update_inits_X.SIS <- function(pars, y0, i) {
  with(pars$ix$X[[i]], {
    S = y0[S_ix]
    I = y0[I_ix]
    pars$Xinits[[i]] = make_Xinits_SIS(pars, list(), S0=S, I0=I)
    return(pars)
})}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SIS <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,I)))
}

#' Plot the density of infected individuals for the SIS model
#'
#' @inheritParams xde_plot_X
#' @export
xde_plot_X.SIS = function(pars, i, clrs=c("darkblue","darkred"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X(vars$XH[[i]], pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SIS model
#'
#' @inheritParams xde_lines_X
#'
#' @export
xde_lines_X.SIS = function(XH, pars, clrs=c("darkblue","darkred"), llty=1){
  with(XH,{
    if(pars$nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 2, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, S[,i], col=clrs[i,1], lty = llty[i])
        lines(time, I[,i], col=clrs[i,2], lty = llty[i])
      }
    }
  })}

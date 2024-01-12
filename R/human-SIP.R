# specialized methods for the human SIP model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIP model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIP <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  X = with(pars$Xpar[[i]], c*I)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIP model.
#' @inheritParams F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIP <- function(t, y, pars, i){
  with(pars$ix$X[[i]],{
    S = y[S_ix]
    I = y[I_ix]
    P = y[P_ix]
    return(S+I+P)
})}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIP model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIP <- function(varslist, pars, i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIP model.
#' @inheritParams F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIP <- function(y, pars,i) {
  with(pars$Xpar[[i]], b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIP model.
#' @inheritParams dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIP <- function(t, y, pars, i){

  foi <- pars$FoI[[i]]

  with(pars$ix$X[[i]], {
    S <- y[S_ix]
    I <- y[I_ix]
    P <- y[P_ix]
    H <- F_H(t, y, pars, i)

   with(pars$Xpar[[i]], {

      dS <- Births(t, H, pars, i)-foi*S -xi*S + r*I + eta*P + dHdt(t, S, pars, i)
      dI <- (1-rho)*foi*S - (r+xi)*I + dHdt(t, I, pars, i)
      dP <- rho*foi*S + xi*(S+I) - eta*P + dHdt(t, P, pars, i)

      return(c(dS, dI, dP))
    })
  })
}


#' @title Compute the HTC for the SIP model
#' @description Implements [HTC] for the SIP model with demography.
#' @inheritParams HTC
#' @return a [numeric] vector
#' @export
HTC.SIP <- function(pars, i) {
  with(pars$Xpar[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}

#' @title Setup Xpar.SIP
#' @description Implements [setup_Xpar] for the SIP model
#' @inheritParams setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SIP = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SIP(pars$nStrata, Xopts)
  return(pars)
}

#' @title Setup Xinits.SIP
#' @description Implements [setup_Xinits] for the SIP model
#' @inheritParams setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIP = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = make_Xinits_SIP(pars$nStrata, Xopts, H0=pars$Hpar[[i]]$H)
  return(pars)
}

#' @title Make parameters for SIP human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param xi background treatment rate
#' @return a [list]
#' @export
make_Xpar_SIP = function(nStrata, Xopts=list(),
                         b=0.55, r=1/180, c=0.15,
                         rho=.1, eta=1/25, xi=1/365){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIP")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$rho = checkIt(rho, nStrata)
    Xpar$eta = checkIt(eta, nStrata)
    Xpar$xi = checkIt(xi, nStrata)

    return(Xpar)
})}


#' @title Make initial values for the SIP human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param H0 the initial human population density
#' @param S0 the initial values of the parameter S
#' @param I0 the initial values of the parameter I
#' @param P0 the initial values of the parameter P
#' @return a [list]
#' @export
make_Xinits_SIP = function(nStrata, Xopts = list(), H0=NULL, S0=NULL,
                           I0=1, P0=0){with(Xopts,{
  if(is.null(S0)) S0=H0-I0-P0
  S = checkIt(S0, nStrata)
  I = checkIt(I0, nStrata)
  P = checkIt(P0, nStrata)
  return(list(S=S, I=I, P=P))
})}

#' @title Parse the output of deSolve and return variables for the SIP model
#' @description Implements [parse_deout_X] for the SIP model
#' @inheritParams parse_deout_X
#' @return none
#' @export
parse_deout_X.SIP <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    I = deout[,I_ix+1]
    P = deout[,P_ix+1]
    H = S+I+P
    return(list(time=time,S=S,I=I,P=P,H=H))
})}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIP model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIP <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(I_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(P_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, P_ix=P_ix)
  return(pars)
})}

#' @title Make parameters for SIP human model
#' @param pars a [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param xi background treatment rate
#' @return none
#' @export
make_parameters_X_SIP <- function(pars, b, c, r, rho, eta, xi){
  stopifnot(is.numeric(b), is.numeric(c), is.numeric(r), is.numeric(rho), is.numeric(eta), is.numeric(xi))
  Xpar <- list()
  class(Xpar) <- c('SIP')
  Xpar$b <- b
  Xpar$c <- c
  Xpar$r <- r
  Xpar$rho <- rho
  Xpar$eta <- eta
  Xpar$xi <- xi
  pars$Xpar[[1]] <- Xpar
  return(pars)
}

#' @title Make inits for SIP human model
#' @param pars a [list]
#' @param S0 size of infected population in each strata
#' @param I0 size of infected population in each strata
#' @param P0 size of population protected by prophylaxis in each strata
#' @return none
#' @export
make_inits_X_SIP <- function(pars, S0, I0, P0) {
  stopifnot(is.numeric(I0), is.numeric(P0), is.numeric(S0))
  pars$Xinits[[1]] = list(S=S0, I=I0, P=P0)
  return(pars)
}

#' @title Update inits for the SIP human model from a vector of states
#' @inheritParams update_inits_X
#' @return none
#' @export
update_inits_X.SIP <- function(pars, y0, i) {
  with(pars,ix$X[[i]],{
  S = y0[S_ix]
  I = y0[I_ix]
  P = y0[P_ix]
  pars$Xinits[[i]] = make_Xinits_SIP(pars$nStrata, list(), S0=S, I0=I, P0=P)
  return(pars)
})}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams get_inits_X
#' @return none
#' @export
get_inits_X.SIP <- function(pars, i){with(pars$Xinits[[i]],{
  c(S, I, P)
})}


#' Plot the density of infected individuals for the SIP model
#'
#' @inheritParams xde_plot_X
#' @export
xde_plot_X.SIP = function(pars, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X_SIP(vars$XH[[i]], pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SIP model
#'
#' @param XH a list with the outputs of parse_deout_X_SIS
#' @param pars a list that defines an `exDE` model (*e.g.*,  generated by `xde_setup()`)
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_SIP = function(XH, pars, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
  with(XH,{
    if(pars$nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, P, col=clrs[3], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 3, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
        lines(time, P[,i], col=clrs[3,i], lty = llty[i])
      }
    }
  })}


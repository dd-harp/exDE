# specialized methods for a basic adult mosquito model

#' @title Set bionomic parameters to baseline
#' @description Implements [MBionomics] for the basicM model
#' @inheritParams MBionomics
#' @return a named [list]
#' @export
MBionomics.basicM <- function(t, y, pars, s) {

  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f <- f0
    pars$MYZpar[[s]]$q <- q0
    pars$MYZpar[[s]]$g <- g0
    pars$MYZpar[[s]]$sigma <- sigma0
    pars$MYZpar[[s]]$nu <- nu0

    return(pars)
})}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the basicM model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.basicM <- function(t, y, pars, s) {
  0*y[pars$ix$MYZ[[s]]$M_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the basic ecology model
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.basicM <- function(t, y, pars, s) {
  M <- y[pars$ix$MYZ[[s]]$M_ix]
  with(pars$MYZpar[[s]], {
    return(M*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the basicM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.basicM <- function(t, y, pars, s){

  Lambda = pars$Lambda[[s]]

  with(pars$ix$MYZ[[s]],{
    M <- y[M_ix]
    P <- y[P_ix]

    with(pars$MYZpar[[s]],{
      Omega <- make_Omega(g, sigma, calK, nPatches)

      dMdt <- Lambda - (Omega %*% M)
      dPdt <- f*(M - P) - (Omega %*% P)

      return(c(dMdt, dPdt))
    })
  })
 }

#' @title Setup MYZpar for the basicM model
#' @description Implements [setup_MYZpar] for the basicM model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.basicM = function(MYZname, pars, s, MYZopts=list(), EIPmod=NULL, calK){
  pars$MYZpar[[s]] = make_MYZpar_basicM(pars$nPatches, MYZopts, calK)
  return(pars)
}

#' @title Make parameters for basicM ODE adult mosquito model
#' @param nPatches the number of patches
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f blood feeding rate
#' @param q human blood feeding fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list] with a configured MYZpar
#' @export
make_MYZpar_basicM = function(nPatches, MYZopts=list(), calK,
                              g=1/12, sigma=1/8,
                              f=0.3, q=0.95,
                              nu=1, eggsPerBatch=60){
  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))
  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "basicM"

    MYZpar$xde <- "ode"
    class(MYZpar$xde) <- "ode"

    MYZpar$nPatches <- nPatches

    MYZpar$g      <- checkIt(g, nPatches)
    MYZpar$sigma  <- checkIt(sigma, nPatches)
    MYZpar$f      <- checkIt(f, nPatches)
    MYZpar$q      <- checkIt(q, nPatches)
    MYZpar$nu     <- checkIt(nu, nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch

    # Store as baseline values
    MYZpar$g0      <- MYZpar$g
    MYZpar$sigma0  <- MYZpar$sigma
    MYZpar$f0      <- MYZpar$f
    MYZpar$q0      <- MYZpar$q
    MYZpar$nu0     <- MYZpar$nu

    MYZpar$calK <- calK

    MYZpar$Omega <- make_Omega(g, sigma, calK, nPatches)

    return(MYZpar)
})}

#' @title Setup the basicM model
#' @description Implements [setup_MYZinits] for the basicM model
#' @inheritParams setup_MYZinits
#' @return a [list] vector
#' @export
setup_MYZinits.basicM = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = make_MYZinits_basicM(pars$nPatches, MYZopts)
  return(pars)
}



#' @title Make inits for basicM adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @return none
#' @export
make_MYZinits_basicM = function(nPatches, MYZopts = list(),
                            M0=5, P0=1){
  with(MYZopts,{
    M = checkIt(M0, nPatches)
    P = checkIt(P0, nPatches)
    return(list(M=M, P=P))
  })}

#' @title Make inits for RM adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return none
#' @export
update_inits_MYZ.basicM <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  P = y0[P_ix]
  pars$MYZinits[[s]] = make_MYZinits_basicM(pars, list(), M0=M, P0=P)
  return(pars)
})}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the basic M model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.basicM <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(M_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(P_ix, 1)

  pars$max_ix = max_ix

  pars$ix$MYZ[[s]] = list(M_ix=M_ix, P_ix=P_ix)

  return(pars)
})}

#' @title Parse the output of deSolve and return variables for the basicM model
#' @description Implements [parse_deout_MYZ] for the basicM model.
#' @inheritParams parse_deout_MYZ
#' @return [list]
#' @export
parse_deout_MYZ.basicM <- function(deout, pars, s) {
  time = deout[,1]
  with(pars$ix$MYZ[[s]],{
    M = deout[,M_ix+1]
    P = deout[,P_ix+1]
    parous = P/M
  return(list(time=time, M=M, P=P, parous=parous))
})}

#' @title Make parameters for a basic adult mosquito model
#' @param pars a [list]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return none
#' @export
make_parameters_MYZ_basicM <- function(pars, g, sigma, f, q, nu, eggsPerBatch, calK) {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  class(MYZpar) <- "basicM"

  MYZpar$g0 <- g
  MYZpar$sigma0 <- sigma
  MYZpar$f0 <- f
  MYZpar$q0 <- q
  MYZpar$nu0 <- nu
  MYZpar$eggsPerBatch <- eggsPerBatch
  MYZpar$calK <- calK

  MYZpar$g <- g
  MYZpar$sigma <- sigma
  MYZpar$f <- f
  MYZpar$q <- q
  MYZpar$nu <- nu

  pars$MYZpar[[1]] <- MYZpar

  return(pars)
}

#' @title Make inits for basicM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_basicM <- function(pars, M0, P0) {
  pars$MYZinits[[1]] = list(M=M0, P=P0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the basicM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.basicM <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(M, P)
})}

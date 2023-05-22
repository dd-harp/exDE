# specialized methods for a basic adult mosquito model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MosquitoBehavior] for the basicM model
#' @inheritParams MosquitoBehavior
#' @return a named [list]
#' @export
MosquitoBehavior.basicM <- function(t, y, pars) {

  pars$MYZpar$f <- pars$MYZpar$f0
  pars$MYZpar$q <- pars$MYZpar$q0
  pars$MYZpar$g <- pars$MYZpar$g0
  pars$MYZpar$sigma <- pars$MYZpar$sigma0
  pars$MYZpar$nu <- pars$MYZpar$nu0

  return(pars)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the basic ecology model
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.basicM <- function(t, y, pars) {
  with(pars$MYZpar, {
    M <- y[M_ix]
    return(M*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the basicM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.basicM <- function(t, y, pars, Lambda, kappa=NULL) {

  nPatches <- pars$nPatches

  with(pars$MYZpar,{

    M <- y[M_ix]
    P <- y[P_ix]

    Omega <- make_Omega(g, sigma, calK, nPatches)

    dMdt <- Lambda - (Omega %*% M)
    dPdt <- f*(M - P) - (Omega %*% P)

    return(c(dMdt, dPdt))
  })
}

#' @title Setup the basicM model for adult mosquitoes
#' @description Implements [setup_MYZ] for the basicM model
#' @inheritParams setup_MYZ
#' @return a [list]
#' @export
setup_MYZ.basicM = function(pars, MYZname,
                           nPatches=1, MYZopts=list(),
                           calK=diag(1)){

  pars$MYZname = "basicM"
  pars$nPatches = nPatches

  pars = make_MYZpar_basicM(pars, MYZopts, calK)
  pars = make_MYZinits_basicM(pars, MYZopts)

  return(pars)
}


#' @title Make parameters for basicM ODE adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param solve_as is either `ode` to solve as an ode or `dde` to solve as a dde
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f blood feeding rate
#' @param q human blood feeding fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @export
make_MYZpar_basicM = function(pars, MYZopts=list(), calK,
                          solve_as = "dde",
                          g=1/12, sigma=1/8,
                          f=0.3, q=0.95,
                          nu=1, eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(pars$nPatches, pars$nPatches))


  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "basicM"

    MYZpar$xde <- "ode"
    class(MYZpar$xde) <- "ode"

    MYZpar$g0      <- checkIt(g, pars$nPatches)
    MYZpar$sigma0  <- checkIt(sigma, pars$nPatches)
    MYZpar$f0      <- checkIt(f, pars$nPatches)
    MYZpar$q0      <- checkIt(q, pars$nPatches)
    MYZpar$nu0     <- checkIt(nu, pars$nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch
    MYZpar$calK <- calK

    pars$MYZpar = MYZpar
    pars = MosquitoBehavior(0, 0, pars)

    return(pars)
  })}

#' @title Make inits for basicM adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @return none
#' @export
make_MYZinits_basicM = function(pars, MYZopts = list(),
                            M0=5, P0=1){
  with(MYZopts,{
    inits = list()
    inits$M0 = checkIt(M0, pars$nPatches)
    inits$P0 = checkIt(P0, pars$nPatches)

    pars$MYZinits = inits
    return(pars)
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the basic M model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.basicM <- function(pars) {

  pars$MYZpar$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$M_ix, 1)

  pars$MYZpar$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$P_ix, 1)

  return(pars)
}

#' @title Parse the output of deSolve and return variables for the basicM model
#' @description Implements [parse_deout_MYZ] for the basicM model.
#' @inheritParams parse_deout_MYZ
#' @return none
#' @export
parse_deout_MYZ.basicM <- function(varslist, deout, pars) {
  varslist$M = deout[,pars$MYZpar$M_ix+1]
  varslist$P = deout[,pars$MYZpar$P_ix+1]
  return(varslist)
}

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

  pars$MYZpar <- MYZpar
  pars = MosquitoBehavior(0, 0, pars)
  return(pars)
}

#' @title Make inits for basicM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_basicM <- function(pars, M0, P0) {
  pars$MYZinits = list(M0=M0, P0=P0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the basicM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.basicM <- function(pars) {with(pars$MYZinits,{
  c(M0, P0)
})}

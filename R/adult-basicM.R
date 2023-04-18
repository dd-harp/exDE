# specialized methods for a basic adult mosquito model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MosquitoBehavior] for the RM model
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
#' @description Implements [dMYZdt] for the RM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.basicM <- function(t, y, pars, Lambda, kappa=NULL) {

  nPatches <- pars$nPatches

  with(pars$MYZpar,{

    M <- y[M_ix]
    P <- y[P_ix]

    Omega <- make_Omega(g, sigma, calK, nPatches)
    Upsilon <- expm(-Omega*eip)

    dMdt <- Lambda - (Omega %*% M)
    dPdt <- f*(M - P) - (Omega %*% P)

    return(c(dMdt, dPdt))
  })
}

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

#' @title Make parameters for a basic adult mosquito model
#' @param pars an [environment]
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

#' @title Make inits for RM adult mosquito model
#' @param pars an [environment]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_basicM <- function(pars, M0, P0) {
  pars$Minits = list(M0=M0, P0=P0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.basicM <- function(pars) {with(pars$Minits,{
  c(M0, P0)
})}

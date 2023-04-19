# specialized methods for the adult mosquito RM model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MosquitoBehavior] for the RM model
#' @inheritParams MosquitoBehavior
#' @return a named [list]
#' @export
MosquitoBehavior.RM <- function(t, y, pars) {

  pars$MYZpar$f <- pars$MYZpar$f0
  pars$MYZpar$q <- pars$MYZpar$q0
  pars$MYZpar$g <- pars$MYZpar$g0
  pars$MYZpar$sigma <- pars$MYZpar$sigma0
  pars$MYZpar$nu <- pars$MYZpar$nu0

  return(pars)
}

#' @title Density of infectious mosquitoes
#' @description Implements [F_Z] for the RM model.
#' @inheritParams F_Z
#' @return a [numeric] vector of length `nPatches`
#' @export
F_Z.RM <- function(t, y, pars) {
  y[pars$MYZpar$Z_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RM <- function(t, y, pars) {
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
dMYZdt.RM_ode <- function(t, y, pars, Lambda, kappa) {

  nPatches <- pars$nPatches

  with(pars$MYZpar,{

    M <- y[M_ix]
    P <- y[P_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]

    Omega <- make_Omega(g, sigma, calK, nPatches)
    Upsilon <- expm(-Omega*eip)

    dMdt <- Lambda - (Omega %*% M)
    dPdt <- f*(M - P) - (Omega %*% P)
    dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
    dZdt <- Upsilon %*% diag(f*q*kappa, nPatches) %*% (M - Y) - (Omega %*% Z)

    return(c(dMdt, dPdt, dYdt, dZdt))
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RM DDE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @importFrom deSolve lagderiv
#' @export
dMYZdt.RM_dde <- function(t, y, pars, Lambda, kappa) {

  nPatches <- pars$nPatches
  eip <- pars$MYZpar$eip

  with(pars$MYZpar,{

    if (t < eip) {
      M_eip <- pars$MYZinits$M0
      Y_eip <- pars$MYZinits$Y0
      kappa_eip <- kappa
      f_eip <- f
      q_eip <- q
      g_eip <- g
      sigma_eip <- sigma
    } else {
      M_eip <- lagvalue(t = t - eip, nr = M_ix)
      Y_eip <- lagvalue(t = t - eip, nr = Y_ix)
      kappa_eip <- lagderiv(t = t-eip, nr = kappa_ix)
      f_eip <- lagderiv(t = t-eip, nr = f_ix)
      q_eip <- lagderiv(t = t-eip, nr = q_ix)
      g_eip <- lagderiv(t = t-eip, nr = g_ix)
      sigma_eip <- lagderiv(t = t-eip, nr = sigma_ix)
    }

    M <- y[M_ix]
    P <- y[P_ix]
    Y <- y[Y_ix]
    Z <- y[Z_ix]
    Upsilon <- matrix(data = y[Upsilon_ix], nrow = nPatches, ncol = nPatches)

    Omega <- make_Omega(g, sigma, calK, nPatches)
    Omega_eip <- make_Omega(g_eip, sigma_eip, calK, nPatches)

    dMdt <- Lambda - (Omega %*% M)
    dPdt <- f*(M - P) - (Omega %*% P)
    dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
    dZdt <- Upsilon %*% diag(f_eip*q_eip*kappa_eip, nPatches) %*% (M_eip - Y_eip) - (Omega %*% Z)
    dUdt <- as.vector((Omega_eip - Omega) %*% Upsilon)

    return(c(dMdt, dPdt, dYdt, dZdt, dUdt, kappa, f, q, g, sigma))
  })
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.RM_ode <- function(pars) {

  pars$MYZpar$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$M_ix, 1)

  pars$MYZpar$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$P_ix, 1)

  pars$MYZpar$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Y_ix, 1)

  pars$MYZpar$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Z_ix, 1)

  return(pars)
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RM model.
#' @inheritParams make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.RM_dde <- function(pars) {

  pars$MYZpar$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$M_ix, 1)

  pars$MYZpar$P_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$P_ix, 1)

  pars$MYZpar$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Y_ix, 1)

  pars$MYZpar$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Z_ix, 1)

  pars$MYZpar$Upsilon_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches^2)
  pars$max_ix <- tail(pars$MYZpar$Upsilon_ix, 1)

  pars$MYZpar$kappa_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$kappa_ix, 1)

  pars$MYZpar$f_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$f_ix, 1)

  pars$MYZpar$q_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$q_ix, 1)

  pars$MYZpar$g_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$g_ix, 1)

  pars$MYZpar$sigma_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$sigma_ix, 1)

  return(pars)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param pars an [environment]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip length of extrinsic incubation period
#' @param solve_as is either `ode` to solve as an ode or `dde` to solve as a dde
#' @return none
#' @export
make_parameters_MYZ_RM <- function(pars, g, sigma, f, q, nu, eggsPerBatch, eip, calK, solve_as = 'dde') {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()

  xde <- solve_as
  class(xde) <- solve_as
  MYZpar$xde <- xde
  if(solve_as == 'dde'){
    class(MYZpar) <- c('RM', 'RM_dde')
    pars$xde <- xde
  }
  else if(solve_as == 'ode') class(MYZpar) <- c('RM', 'RM_ode')

  MYZpar$g0 <- g
  MYZpar$sigma0 <- sigma
  MYZpar$f0 <- f
  MYZpar$q0 <- q
  MYZpar$nu0 <- nu
  MYZpar$eggsPerBatch <- eggsPerBatch
  MYZpar$eip <- eip
  MYZpar$calK <- calK

  pars$MYZpar <- MYZpar
  pars = MosquitoBehavior.RM(0, 0, pars)
  return(pars)
}

#' @title Make inits for RM adult mosquito model
#' @param pars an [environment]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_RM_ode <- function(pars, M0, P0, Y0, Z0) {
  pars$MYZinits = list(M0=M0, P0=P0, Y0=Y0, Z0=Z0)
  return(pars)
}

#' @title Make inits for RM adult mosquito model
#' @param pars an [environment]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @param Upsilon0 the initial values of Upsilon
#' @return none
#' @export
make_inits_MYZ_RM_dde <- function(pars, M0, P0, Y0, Z0, Upsilon0) {
  pars$MYZinits = list(M0=M0, P0=P0, Y0=Y0, Z0=Z0, Upsilon0=Upsilon0, rep(0, 5*pars$nPatches))
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.RM_ode <- function(pars) {with(pars$MYZinits,{
  c(M0, P0, Y0, Z0)
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.RM_dde <- function(pars) {with(pars$MYZinits,{
  c(M0, P0, Y0, Z0, as.vector(Upsilon0), rep(0, 5*pars$nPatches))
})}

# specialized methods for the adult mosquito RM model

#' @title Entomological inoculation rate on human strata
#' @description Implements [F_EIR] for the generalized RM model.
#' @inheritParams F_EIR
#' @return a [numeric] vector of length `nStrata`
#' @export
F_EIR.RM <- function(t, y, pars) {
  Z <- y[pars$Z_ix]
  with(pars$MYZpar, {
    return(
      as.vector(pars$beta %*% diag(f*q, nrow = pars$nPatches) %*% Z)
    )
  })
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_kappa] for the generalized RM ODE model.
#' @inheritParams F_kappa
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa.RM_ode <- function(t, y, pars) {
  x <- F_x(t, y, pars)
  as.vector(pars$betaT %*% x)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_kappa] for the generalized RM DDE model.
#' @inheritParams F_kappa
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa.RM_dde <- function(t, y, pars) {

  x <- F_x(t, y, pars)
  x_tau <- F_x_tau(t, y, pars, pars$MYZpar$tau)

  kappa <- matrix(data = 0, nrow = 2, ncol = pars$nPatches)
  kappa[1, ] <- as.vector(pars$betaT %*% x)
  kappa[2, ] <- as.vector(pars$betaT %*% x_tau)
  return(kappa)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the generalized RM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RM <- function(t, y, pars) {
  G <- y[pars$G_ix]
  with(pars$MYZpar, {
    return(G*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the generalized RM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RM_ode <- function(t, y, pars, Lambda, kappa) {

  M <- y[pars$M_ix]
  G <- y[pars$G_ix]
  Y <- y[pars$Y_ix]
  Z <- y[pars$Z_ix]

  with(pars$MYZpar, {
    dMdt <- Lambda - (Omega %*% M)
    dGdt <- diag(f, pars$nPatches, pars$nPatches) %*% (M - G) - (nu * G) - (Omega %*% G)
    dYdt <- diag(f*q*kappa) %*% (M - Y) - (Omega %*% Y)
    dZdt <- OmegaEIP %*% diag(f*q*kappa) %*% (M - Y) - (Omega %*% Z)
    return(c(dMdt, dGdt, dYdt, dZdt))
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the generalized RM DDE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @export
dMYZdt.RM_dde <- function(t, y, pars, Lambda, kappa) {
  kappa_t <- kappa[1, ]
  kappa_tau <- kappa[2, ]

  M <- y[pars$M_ix]
  G <- y[pars$G_ix]
  Y <- y[pars$Y_ix]
  Z <- y[pars$Z_ix]

  tau <- pars$MYZpar$tau

  if (t < tau) {
    M_tau <- pars$MYZpar$M0
    Y_tau <- pars$MYZpar$Y0
  } else {
    M_tau <- lagvalue(t = t - tau, nr = pars$M_ix)
    Y_tau <- lagvalue(t = t - tau, nr = pars$Y_ix)
  }

  with(pars$MYZpar, {
    dMdt <- Lambda - (Omega %*% M)
    dGdt <- diag(f, pars$nPatches, pars$nPatches) %*% (M - G) - (nu * G) - (Omega %*% G)
    dYdt <- diag(f*q*kappa_t) %*% (M - Y) - (Omega %*% Y)
    dZdt <- OmegaEIP %*% diag(f*q*kappa_tau) %*% (M_tau - Y_tau) - (Omega %*% Z)
    return(c(dMdt, dGdt, dYdt, dZdt))
  })
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_index_MYZ] for the generalized RM model.
#' @inheritParams make_index_MYZ
#' @return the modified parameter [list]
#' @importFrom utils tail
#' @export
make_index_MYZ.RM <- function(pars) {
  pars$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$M_ix, 1)

  pars$G_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$G_ix, 1)

  pars$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$Y_ix, 1)

  pars$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$Z_ix, 1)
  return(pars)
}

#' @noRd
make_parameters_MYZ_RM <- function(MYZpar, Omega, OmegaEIP, f, q, nu, eggsPerBatch, M0, G0, Y0, Z0) {
  stopifnot(vapply(as.list(match.call()), is.numeric, logical(1)))
  stopifnot(inherits(Omega, 'matrix'), inherits(OmegaEIP, 'matrix'))
  # stopifnot(is.numeric(psi), is.numeric(phi), is.numeric(theta), is.numeric(L0))
  MYZpar$Omega <- Omega
  MYZpar$OmegaEIP <- OmegaEIP
  MYZpar$f <- f
  MYZpar$q <- q
  MYZpar$nu <- nu
  MYZpar$eggsPerBatch <- eggsPerBatch
  MYZpar$M0 <- M0
  MYZpar$G0 <- G0
  MYZpar$Y0 <- Y0
  MYZpar$Z0 <- Z0
  return(MYZpar)
}

#' @title Make parameters for generalized RM ODE adult mosquito model
#' @param Omega mosquito demography matrix
#' @param OmegaEIP mosquito demography matrix through the EIP
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate of gravid mosquitoes
#' @param eggsPerBatch eggs laid per oviposition
#' @param M0 total mosquito density at each patch
#' @param G0 gravid mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list] with classes `RM`, `RM_ode`.
#' @export
make_parameters_MYZ_RM_ode <- function(Omega, OmegaEIP, f, q, nu, eggsPerBatch, M0, G0, Y0, Z0) {
  MYZpar <- list()
  class(MYZpar) <- c('RM', 'RM_ode')
  MYZpar <- make_parameters_MYZ_RM(MYZpar = MYZpar, Omega = Omega, OmegaEIP = OmegaEIP, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = M0, G0 = G0, Y0 = Y0, Z0 = Z0)
  return(MYZpar)
}

#' @title Make parameters for generalized RM DDE adult mosquito model
#' @param Omega mosquito demography matrix
#' @param OmegaEIP mosquito demography matrix through the EIP
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate of gravid mosquitoes
#' @param eggsPerBatch eggs laid per oviposition
#' @param M0 total mosquito density at each patch
#' @param G0 gravid mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @param tau length of extrinsic incubation period
#' @return a [list] with classes `RM`, `RM_dde`.
#' @export
make_parameters_MYZ_RM_dde <- function(Omega, OmegaEIP, f, q, nu, eggsPerBatch, tau, M0, G0, Y0, Z0) {
  MYZpar <- list()
  class(MYZpar) <- c('RM', 'RM_dde')
  MYZpar <- make_parameters_MYZ_RM(MYZpar = MYZpar, Omega = Omega, OmegaEIP = OmegaEIP, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = M0, G0 = G0, Y0 = Y0, Z0 = Z0)
  MYZpar$tau <- tau
  return(MYZpar)
}

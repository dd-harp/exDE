# specialized methods for the adult mosquito GeRM model

#' @title Compute bloodfeeding and mortality rates
#' @description Implements [MosquitoBehavior] for the generalized GeRM model.
#' @inheritParams MosquitoBehavior
#' @return a named [list]
#' @export
MosquitoBehavior.GeRM <- function(t, y, pars) {

  MosyBehavior <- list()
  MosyBehavior$f <- rep(pars$MYZpar$f, 2)
  attr(MosyBehavior$f, 'time') <- c(t, t - pars$MYZpar$tau)
  MosyBehavior$q <- rep(pars$MYZpar$q, 2)
  MosyBehavior$g <- rep(pars$MYZpar$g, 2)

  return(MosyBehavior)
}

#' @title Time spent host seeking/feeding and resting/ovipositing
#' @description Implements [F_tau] for the generalized GeRM model.
#' @inheritParams F_tau
#' @return [NULL]
#' @export
F_tau.GeRM <- function(t, y, pars) {
  NULL
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_kappa] for the generalized GeRM ODE model.
#' @inheritParams F_kappa
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa.GeRM_ode <- function(t, y, pars) {
  x <- F_x(t, y, pars)
  beta <- F_beta(t, y, pars)
  as.vector(t(beta) %*% x)
}

#' @title Net infectiousness of human population to mosquitoes
#' @description Implements [F_kappa] for the generalized GeRM DDE model.
#' @inheritParams F_kappa
#' @return a [numeric] vector of length `nPatches`
#' @export
F_kappa.GeRM_dde <- function(t, y, pars) {

  x <- F_x(t, y, pars)
  x_lag <- F_x_lag(t, y, pars, pars$MYZpar$tau)

  beta <- F_beta(t, y, pars)
  beta_lag <- F_beta_lag(t, y, pars, pars$MYZpar$tau)

  kappa <- matrix(data = 0, nrow = 2, ncol = pars$nPatches)
  kappa[1, ] <- as.vector(t(beta) %*% x)
  kappa[2, ] <- as.vector(t(beta_lag) %*% x_lag)
  return(kappa)
}

#' @title Density of infectious mosquitoes
#' @description Implements [F_Z] for the generalized GeRM model.
#' @inheritParams F_Z
#' @return a [numeric] vector of length `nPatches`
#' @export
F_Z.GeRM <- function(t, y, pars) {
  y[pars$Z_ix]
}

#' @title Density of lagged infectious mosquitoes
#' @description Implements [F_Z_lag] for the generalized GeRM model.
#' @inheritParams F_Z_lag
#' @return a [numeric] vector of length `nPatches`
#' @importFrom deSolve lagvalue
#' @export
F_Z_lag.GeRM <- function(t, y, pars, lag) {
  if (t < lag) {
    return(pars$MYZpar$Z0)
  } else {
    return(lagvalue(t = t - lag, nr = pars$Z_ix))
  }
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the generalized GeRM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.GeRM <- function(t, y, pars) {
  G <- y[pars$G_ix]
  with(pars$MYZpar, {
    return(G*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the generalized GeRM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.GeRM_ode <- function(t, y, pars, Lambda, kappa, MosyBehavior) {

  nPatches <- pars$nPatches

  M <- y[pars$M_ix]
  G <- y[pars$G_ix]
  Y <- y[pars$Y_ix]
  Z <- y[pars$Z_ix]
  Upsilon <- matrix(data = y[pars$Upsilon_ix], nrow = nPatches, ncol = nPatches)

  f <- MosyBehavior$f
  q <- MosyBehavior$q
  g <- MosyBehavior$g

  Omega <- make_Omega(g = g[1], sigma = pars$MYZpar$sigma, K = pars$MYZpar$calK, nPatches = nPatches)
  Omega_eip <- make_Omega(g = g[2], sigma = pars$MYZpar$sigma, K = pars$MYZpar$calK, nPatches = nPatches)

  dMdt <- Lambda - (Omega %*% M)
  dGdt <- diag(f[1], nPatches) %*% (M - G) - (pars$MYZpar$nu * G) - (Omega %*% G)
  dYdt <- diag(f[1]*q[1]*kappa, nPatches) %*% (M - Y) - (Omega %*% Y)
  dZdt <- Upsilon %*% diag(f[2]*q[2]*kappa, nPatches) %*% (M - Y) - (Omega %*% Z)
  dUdt <- as.vector((Omega_eip - Omega) %*% Upsilon)
  return(c(dMdt, dGdt, dYdt, dZdt, dUdt))
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the generalized GeRM DDE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @export
dMYZdt.GeRM_dde <- function(t, y, pars, Lambda, kappa, MosyBehavior) {
  kappa_t <- kappa[1, ]
  kappa_eip <- kappa[2, ]

  nPatches <- pars$nPatches

  M <- y[pars$M_ix]
  G <- y[pars$G_ix]
  Y <- y[pars$Y_ix]
  Z <- y[pars$Z_ix]
  Upsilon <- matrix(data = y[pars$Upsilon_ix], nrow = nPatches, ncol = nPatches)

  tau <- pars$MYZpar$tau

  if (t < tau) {
    M_tau <- pars$MYZpar$M0
    Y_tau <- pars$MYZpar$Y0
  } else {
    M_tau <- lagvalue(t = t - tau, nr = pars$M_ix)
    Y_tau <- lagvalue(t = t - tau, nr = pars$Y_ix)
  }

  f <- MosyBehavior$f
  q <- MosyBehavior$q
  g <- MosyBehavior$g

  Omega <- make_Omega(g = g[1], sigma = pars$MYZpar$sigma, K = pars$MYZpar$calK, nPatches = nPatches)
  Omega_eip <- make_Omega(g = g[2], sigma = pars$MYZpar$sigma, K = pars$MYZpar$calK, nPatches = nPatches)

  dMdt <- Lambda - (Omega %*% M)
  dGdt <- diag(f[1], nPatches) %*% (M - G) - (pars$MYZpar$nu * G) - (Omega %*% G)
  dYdt <- diag(f[1]*q[1]*kappa_t, nPatches) %*% (M - Y) - (Omega %*% Y)
  dZdt <- Upsilon %*% diag(f[2]*q[2]*kappa_eip, nPatches) %*% (M_tau - Y_tau) - (Omega %*% Z)
  dUdt <- as.vector((Omega_eip - Omega) %*% Upsilon)
  return(c(dMdt, dGdt, dYdt, dZdt, dUdt))
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_index_MYZ] for the generalized GeRM model.
#' @inheritParams make_index_MYZ
#' @return the modified parameter [list]
#' @importFrom utils tail
#' @export
make_index_MYZ.GeRM <- function(pars) {
  pars$M_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$M_ix, 1)

  pars$G_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$G_ix, 1)

  pars$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$Y_ix, 1)

  pars$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$Z_ix, 1)

  pars$Upsilon_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches^2)
  pars$max_ix <- tail(pars$Upsilon_ix, 1)
  return(pars)
}


#' @noRd
make_parameters_MYZ_GeRM <- function(MYZpar, g, sigma, calK, f, q, nu, eggsPerBatch, tau, M0, G0, Y0, Z0) {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch), is.numeric(M0), is.numeric(G0), is.numeric(Y0), is.numeric(Z0))
  MYZpar$g <- g
  MYZpar$sigma <- sigma
  MYZpar$calK <- calK
  MYZpar$f <- f
  MYZpar$q <- q
  MYZpar$nu <- nu
  MYZpar$eggsPerBatch <- eggsPerBatch
  MYZpar$tau <- tau
  MYZpar$M0 <- M0
  MYZpar$G0 <- G0
  MYZpar$Y0 <- Y0
  MYZpar$Z0 <- Z0
  return(MYZpar)
}

#' @title Make parameters for generalized GeRM ODE adult mosquito model
#' @param pars an [environment]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate of gravid mosquitoes
#' @param eggsPerBatch eggs laid per oviposition
#' @param tau length of extrinsic incubation period
#' @param M0 total mosquito density at each patch
#' @param G0 gravid mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return none
#' @export
make_parameters_MYZ_GeRM_ode <- function(pars, g, sigma, calK, f, q, nu, eggsPerBatch, tau, M0, G0, Y0, Z0) {
  stopifnot(is.environment(pars))
  stopifnot(nrow(calK) == pars$nPatches && ncol(calK) == pars$nPatches)
  MYZpar <- list()
  class(MYZpar) <- c('GeRM', 'GeRM_ode')
  MYZpar <- make_parameters_MYZ_GeRM(MYZpar = MYZpar, g = g, sigma = sigma, calK = calK, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, tau = tau, M0 = M0, G0 = G0, Y0 = Y0, Z0 = Z0)
  pars$MYZpar <- MYZpar
}

#' @title Make parameters for generalized GeRM DDE adult mosquito model
#' @param pars an [environment]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate of gravid mosquitoes
#' @param eggsPerBatch eggs laid per oviposition
#' @param tau length of extrinsic incubation period
#' @param M0 total mosquito density at each patch
#' @param G0 gravid mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return none
#' @export
make_parameters_MYZ_GeRM_dde <- function(pars, g, sigma, calK, f, q, nu, eggsPerBatch, tau, M0, G0, Y0, Z0) {
  MYZpar <- list()
  class(MYZpar) <- c('GeRM', 'GeRM_dde')
  MYZpar <- make_parameters_MYZ_GeRM(MYZpar = MYZpar, g = g, sigma = sigma, calK = calK, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, tau = tau, M0 = M0, G0 = G0, Y0 = Y0, Z0 = Z0)
  pars$MYZpar <- MYZpar
}

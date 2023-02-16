library(expm)
library(MASS)
library(deSolve)

numeric_tol <- 1e-5

test_that("test equilibrium with RM adults (ODE), basic competition", {

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  wf <- rep(1, nStrata)

  f <- 0.3
  q <- 0.9
  g <- 1/10
  sigma <- 1/100
  nu <- 1/2
  eggsPerBatch <- 30

  tau <- 11

  # mosquito movement calK
  calK <- matrix(0, nPatches, nPatches)
  calK[upper.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK[lower.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  # omega matrix
  Omega <- make_Omega(g, sigma, calK, nPatches)
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * tau)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  # human PfPR and H
  pfpr <- rep(0.3, times = nStrata)
  H <- rpois(n = nStrata, lambda = 1000)
  X <- rbinom(n = nStrata, size = H, prob = pfpr)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nPatches, nPatches) %*% ((r*X) / (H - X))

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  kappa <- t(beta) %*% (X*c)

  # equilibrium solutions for adults
  Z <- diag(1/(f*q), nPatches, nPatches) %*% ginv(beta) %*% EIR
  MY <- diag(1/as.vector(f*q*kappa), nPatches, nPatches) %*% OmegaEIP_inv %*% Omega %*% Z
  Y <- Omega_inv %*% (diag(as.vector(f*q*kappa), nPatches, nPatches) %*% MY)
  M <- MY + Y
  G <- solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M

  # equilibrium solutions for aquatic
  calN <- matrix(0, nPatches, nHabitats)
  diag(calN) <- 1

  calU <- matrix(0, nHabitats, nPatches)
  diag(calU) <- 1

  alpha <- as.vector(solve(calN) %*% Lambda)

  psi <- 1/10
  phi <- 1/12
  eta <- as.vector(calU %*% G * nu * eggsPerBatch)

  L <- alpha/psi
  theta <- (eta - psi*L - phi*L)/(L^2)

  # parameters for exDE
  params <- list()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN

  params = make_parameters_MYZ_GeRM(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as = "ode")
  params = make_inits_MYZ_GeRM(pars = params, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z), Upsilon0=OmegaEIP)
  params = make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta = theta)
  params = make_inits_L_basic(pars = params,  L0 = L)
  params = make_parameters_vc_null(pars = params)
  params = make_parameters_exogenous_null(pars = params)

  params = make_indices(params)

  # set initial conditions
  y0 <- get_inits(params)

  # mimic MosyBehavior
  MosyBehavior <- list()
  MosyBehavior$f <- rep(params$MYZpar$f, 2)
  attr(MosyBehavior$f, 'time') <- c(0, 0 - params$MYZpar$tau)
  MosyBehavior$q <- rep(params$MYZpar$q, 2)
  MosyBehavior$g <- rep(params$MYZpar$g, 2)

  # run simulation
  out <- deSolve::ode(y = y0, times = c(0,50), func = xDE_diffeqn_mosy, parms = params, method = "lsoda", kappa = as.vector(kappa), MosyBehavior = MosyBehavior)

  expect_equal(as.vector(out[2, params$L_ix+1]), as.vector(L), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
})

test_that("test equilibrium with RM adults (DDE), basic competition", {

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  wf <- rep(1, nStrata)

  f <- 0.3
  q <- 0.9
  g <- 1/10
  sigma <- 1/100
  nu <- 1/2
  eggsPerBatch <- 30

  tau <- 11

  # mosquito movement calK
  calK <- matrix(0, nPatches, nPatches)
  calK[upper.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK[lower.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  # omega matrix
  Omega <- make_Omega(g, sigma, calK, nPatches)
  Omega_inv <- solve(Omega)
  OmegaEIP <- expm::expm(-Omega * tau)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  # human PfPR and H
  pfpr <- rep(0.3, times = nStrata)
  H <- rpois(n = nStrata, lambda = 1000)
  X <- rbinom(n = nStrata, size = H, prob = pfpr)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nPatches, nPatches) %*% ((r*X) / (H - X))

  # ambient pop
  W <- Psi %*% H

  # biting distribution matrix
  beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), nPatches , nPatches)

  # kappa
  kappa <- t(beta) %*% (X*c)

  # equilibrium solutions for adults
  Z <- diag(1/(f*q), nPatches, nPatches) %*% ginv(beta) %*% EIR
  MY <- diag(1/as.vector(f*q*kappa), nPatches, nPatches) %*% OmegaEIP_inv %*% Omega %*% Z
  Y <- Omega_inv %*% (diag(as.vector(f*q*kappa), nPatches, nPatches) %*% MY)
  M <- MY + Y
  G <- solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M

  # equilibrium solutions for aquatic
  calN <- matrix(0, nPatches, nHabitats)
  diag(calN) <- 1

  calU <- matrix(0, nHabitats, nPatches)
  diag(calU) <- 1

  alpha <- as.vector(solve(calN) %*% Lambda)

  psi <- 1/10
  phi <- 1/12
  eta <- as.vector(calU %*% G * nu * eggsPerBatch)

  L <- alpha/psi
  theta <- (eta - psi*L - phi*L)/(L^2)

  params <- list()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN

  # parameters for exDE
  params = make_parameters_MYZ_GeRM(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch)
  params = make_inits_MYZ_GeRM(pars = params, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z), Upsilon0=OmegaEIP)
  params = make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta = theta)
  params = make_inits_L_basic(pars = params,  L0 = L)
  params = make_parameters_vc_null(pars = params)
  params = make_parameters_exogenous_null(pars = params)

  params = make_indices(params)

  # set initial conditions
  y0 <- get_inits(params)

  # mimic MosyBehavior
  MosyBehavior <- list()
  MosyBehavior$f <- rep(params$MYZpar$f, 2)
  attr(MosyBehavior$f, 'time') <- c(0, 0 - params$MYZpar$tau)
  MosyBehavior$q <- rep(params$MYZpar$q, 2)
  MosyBehavior$g <- rep(params$MYZpar$g, 2)

  # run simulation
  out <- deSolve::dede(y = y0, times = c(0,350), func = xDE_diffeqn_mosy, parms = params, method = "lsoda", kappa = t(cbind(kappa,kappa)), MosyBehavior = MosyBehavior)

  expect_equal(as.vector(out[2, params$L_ix+1]), as.vector(L), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
})

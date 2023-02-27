library(expm)
library(MASS)
library(deSolve)

numeric_tol <- 1e-5

test_that("test equilibrium with RM adults (ODE), SIP humans, trace", {

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  eta <- c(1/30, 1/40)
  rho <- c(0.05, 0.1)
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
  membershipH = c(1,2)
  searchWtsH = c(1,1)
  X <- rbinom(n = nStrata, size = H, prob = pfpr)
  P <- diag(1/eta) %*% diag(rho/(1-rho)) %*% (r*X)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nStrata) %*% diag(1/(1-rho)) %*% ((r*X)/(H-X-P))

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

  # parameters for exDE
  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN

  params = make_parameters_MYZ_GeRM(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as="ode")
  params = make_inits_MYZ_GeRM(pars = params, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z), Upsilon0=OmegaEIP)
  params = make_parameters_demography_null(pars = params, H=H, membershipH=membershipH,
                                           searchWtsH=searchWtsH, TimeSpent=Psi)
  params = make_parameters_X_SIP(pars = params, b = b, c = c, r = r, eta=eta, rho=rho)
  params = make_inits_X_SIP(pars = params, X, P)
  params = make_parameters_L_trace(pars = params, Lambda = as.vector(Lambda))
  params = make_parameters_exogenous_null(pars = params)
  params = make_parameters_vc_null(pars = params)

  params = make_indices(params)


  # set initial conditions
  y0 <- get_inits(params)

  # run simulation
  out <- deSolve::ode(y = y0, times = c(0,50), func = xDE_diffeqn, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$X_ix+1]), as.vector(X), tolerance = numeric_tol)
})

test_that("test equilibrium with RM adults (DDE), SIP humans, trace", {

  # set number of patches and strata
  nPatches <- 2
  nStrata <- 2
  nHabitats <- 2

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  eta <- c(1/30, 1/40)
  rho <- c(0.05, 0.1)
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
  membershipH = c(1,2)
  searchWtsH = c(1,1)
  X <- rbinom(n = nStrata, size = H, prob = pfpr)
  P <- diag(1/eta) %*% diag(rho/(1-rho)) %*% (r*X)

  # TaR
  Psi <- matrix(rexp(n = nStrata*nPatches), nStrata, nPatches)
  Psi <- Psi/rowSums(Psi)
  Psi <- t(Psi)

  # derived EIR to sustain equilibrium pfpr
  EIR <- diag(1/b, nStrata) %*% diag(1/(1-rho)) %*% ((r*X)/(H-X-P))

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

  # parameters for exDE
  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- nPatches
  params$nHabitats <- nHabitats
  params$calU <- calU
  params$calN <- calN

  params = make_parameters_MYZ_GeRM(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as="ode")
  params = make_inits_MYZ_GeRM(pars = params, M0 = as.vector(M), G0 = as.vector(G), Y0 = as.vector(Y), Z0 = as.vector(Z), Upsilon0=OmegaEIP)
  params = make_parameters_demography_null(pars = params, H=H, membershipH=membershipH,
                                           searchWtsH=searchWtsH, TimeSpent=Psi)
  params = make_parameters_X_SIP(pars = params, b = b, c = c, r = r, eta=eta, rho=rho)
  params = make_inits_X_SIP(pars = params, X, P)
  params = make_parameters_L_trace(pars = params, Lambda = as.vector(Lambda))
  params = make_parameters_exogenous_null(pars = params)
  params = make_parameters_vc_null(pars = params)

  params = make_indices(params)


  # set initial conditions
  y0 <- get_inits(params)

  # run simulation
  out <- deSolve::dede(y = y0, times = c(0,50), func = xDE_diffeqn, parms = params, method = "lsoda")

  expect_equal(as.vector(out[2, params$M_ix+1]), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$G_ix+1]), as.vector(G), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Y_ix+1]), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$Z_ix+1]), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(out[2, params$X_ix+1]), as.vector(X), tolerance = numeric_tol)
})

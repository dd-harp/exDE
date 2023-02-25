library(expm)
library(deSolve)


test_that("RM models reach equilibrium", {

  # tolerance for tests comparing floats
  numeric_tol <- 1e-4

  nPatches <- 3
  f <- 0.3
  q <- 0.9
  g <- 1/20
  sigma <- 1/10
  tau <- 11
  nu <- 1/2
  eggsPerBatch <- 30

  calK <- matrix(0, nPatches, nPatches)
  calK[1, 2:3] <- c(0.2, 0.8)
  calK[2, c(1,3)] <- c(0.5, 0.5)
  calK[3, 1:2] <- c(0.7, 0.3)
  calK <- t(calK)

  Omega <- make_Omega(g, sigma, calK, nPatches)
  OmegaEIP <- expm::expm(-Omega * tau)

  kappa <- c(0.1, 0.075, 0.025)
  Lambda <- c(5, 10, 8)

  params <-make_parameters_xde()
  params$nPatches = nPatches

  # ODE
  params = make_parameters_MYZ_GeRM(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as = "ode")
  params = make_inits_MYZ_GeRM(pars = params, M0 = rep(0, nPatches), G0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 =rep(0, nPatches), Upsilon0=OmegaEIP)

  params = make_indices(params)

  # mimic MosyBehavior
  MosyBehavior <- list()
  MosyBehavior$f <- rep(params$MYZpar$f, 2)
  attr(MosyBehavior$f, 'time') <- c(0, 0 - params$MYZpar$tau)
  MosyBehavior$q <- rep(params$MYZpar$q, 2)
  MosyBehavior$g <- rep(params$MYZpar$g, 2)

  # make indices and set up initial conditions
  y0 <- get_inits(params)

  # solve ODEs
  out <- deSolve::ode(y = y0, times = c(0, 730), func = function(t, y, pars, Lambda, kappa, MosyBehavior) {
    list(dMYZdt(t, y, pars, Lambda, kappa, MosyBehavior))
  }, parms = params, method = 'lsoda', Lambda = Lambda, kappa = kappa, MosyBehavior = MosyBehavior)

  # equilibrium solutions (forward)
  Omega_inv <- solve(Omega)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  M_eq <- as.vector(Omega_inv %*% Lambda)
  M_sim <- as.vector(out[2, params$M_ix+1])

  G_eq <- as.vector(solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  G_sim <- as.vector(out[2, params$G_ix+1])

  Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Y_sim <- as.vector(out[2, params$Y_ix+1])

  Z_eq <- as.vector(Omega_inv %*% OmegaEIP %*% diag(f*q*kappa) %*% (M_eq - Y_eq))
  Z_sim <- as.vector(out[2, params$Z_ix+1])

  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(G_eq, G_sim, tolerance = numeric_tol)
  expect_equal(Y_eq, Y_sim, tolerance = numeric_tol)
  expect_equal(Z_eq, Z_sim, tolerance = numeric_tol)

  # equilibrium solutions (backward)
  MY_eq <- as.vector(diag(1/(f*q*kappa)) %*% OmegaEIP_inv %*% Omega %*% Z_eq)
  MY_sim <- M_sim - Y_sim

  Y_eq <- as.vector(Omega_inv %*% diag(f*q*kappa) %*% MY_eq)

  M_eq <- MY_eq + Y_eq

  Lambda_eq <- as.vector(Omega %*% M_eq)

  expect_equal(MY_eq, MY_sim, tolerance = numeric_tol)
  expect_equal(Y_eq, Y_sim, tolerance = numeric_tol)
  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(Lambda_eq, Lambda, tolerance = numeric_tol)

  # DDE
  params = make_parameters_MYZ_GeRM(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch)
  params = make_inits_MYZ_GeRM(pars = params, M0 = rep(0, nPatches), G0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 =rep(0, nPatches), Upsilon0=OmegaEIP)
  params = make_indices(params)

  # solve DDEs
  out <- deSolve::dede(y = y0, times = c(0, 365), func = function(t, y, pars, Lambda, kappa, MosyBehavior) {
    list(dMYZdt(t, y, pars, Lambda, kappa, MosyBehavior))
  }, parms = params, method = 'lsoda', Lambda = Lambda, kappa = rbind(kappa, kappa), MosyBehavior = MosyBehavior)

  # equilibrium solutions (forward)
  M_eq <- as.vector(Omega_inv %*% Lambda)
  M_sim <- as.vector(out[2, params$M_ix+1])

  G_eq <- as.vector(solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  G_sim <- as.vector(out[2, params$G_ix+1])

  Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Y_sim <- as.vector(out[2, params$Y_ix+1])

  Z_eq <- as.vector(Omega_inv %*% OmegaEIP %*% diag(f*q*kappa) %*% (M_eq - Y_eq))
  Z_sim <- as.vector(out[2, params$Z_ix+1])

  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(G_eq, G_sim, tolerance = numeric_tol)
  expect_equal(Y_eq, Y_sim, tolerance = numeric_tol)
  expect_equal(Z_eq, Z_sim, tolerance = numeric_tol)

  # equilibrium solutions (backward)
  MY_eq <- as.vector(diag(1/(f*q*kappa)) %*% OmegaEIP_inv %*% Omega %*% Z_eq)
  MY_sim <- M_sim - Y_sim

  Y_eq <- as.vector(Omega_inv %*% diag(f*q*kappa) %*% MY_eq)

  M_eq <- MY_eq + Y_eq

  Lambda_eq <- as.vector(Omega %*% M_eq)

  expect_equal(MY_eq, MY_sim, tolerance = numeric_tol)
  expect_equal(Y_eq, Y_sim, tolerance = numeric_tol)
  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(Lambda_eq, Lambda, tolerance = numeric_tol)
})

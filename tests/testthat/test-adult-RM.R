library(expm)
library(deSolve)


test_that("RM models reach equilibrium", {

  # tolerance for tests comparing floats
  numeric_tol <- 1e-4

  nPatches <- 3
  f <- rep(0.3, nPatches)
  q <- rep(0.9, nPatches)
  g <- rep(1/20, nPatches)
  sigma <- rep(1/10, nPatches)
  eip <- 11
  nu <- 1/2
  eggsPerBatch <- 30

  calK <- matrix(0, nPatches, nPatches)
  calK[1, 2:3] <- c(0.2, 0.8)
  calK[2, c(1,3)] <- c(0.5, 0.5)
  calK[3, 1:2] <- c(0.7, 0.3)
  calK <- t(calK)

  Omega <- make_Omega(g, sigma, calK, nPatches)
  OmegaEIP <- expm::expm(-Omega * eip)

  kappa <- c(0.1, 0.075, 0.025)
  Lambda <- c(5, 10, 8)

  params <-make_parameters_xde()
  params$nPatches = nPatches

  # ODE
  params = make_parameters_MYZ_RM(pars = params, g = g, sigma = sigma, calK = calK, eip = eip, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as = "ode")
  params = make_inits_MYZ_RM_ode(pars = params, M0 = rep(0, nPatches), P0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 =rep(0, nPatches))

  params = make_indices(params)

  # make indices and set up initial conditions
  y0 <- get_inits(params)

  # solve ODEs
  out <- deSolve::ode(y = y0, times = c(0, 730), func = function(t, y, pars, Lambda, kappa) {
    list(dMYZdt(t, y, pars, Lambda, kappa))
  }, parms = params, method = 'lsoda', Lambda = Lambda, kappa = kappa)

  # equilibrium solutions (forward)
  Omega_inv <- solve(Omega)
  OmegaEIP_inv <- expm::expm(Omega * eip)

  M_eq <- as.vector(Omega_inv %*% Lambda)
  M_sim <- as.vector(out[2, params$ix$MYZ$M_ix+1])

  P_eq <- as.vector(solve(diag(f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  P_sim <- as.vector(out[2, params$ix$MYZ$P_ix+1])

  Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Y_sim <- as.vector(out[2, params$ix$MYZ$Y_ix+1])

  Z_eq <- as.vector(Omega_inv %*% OmegaEIP %*% diag(f*q*kappa) %*% (M_eq - Y_eq))
  Z_sim <- as.vector(out[2, params$ix$MYZ$Z_ix+1])

  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(P_eq, P_sim, tolerance = numeric_tol)
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
  expect_equal(P_eq, P_sim, tolerance = numeric_tol)
  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(Lambda_eq, Lambda, tolerance = numeric_tol)

  # DDE
  params = make_parameters_MYZ_RM(pars = params, g = g, sigma = sigma, calK = calK, eip = eip, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, solve_as = "dde")
  params = make_inits_MYZ_RM_dde(pars = params, M0 = rep(0, nPatches), P0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 =rep(0, nPatches), Upsilon0=as.vector(OmegaEIP))

  params = make_indices(params)

  y0 <- get_inits(params)

  # solve DDEs
  out <- deSolve::dede(y = y0, times = c(0, 365), func = function(t, y, pars, Lambda, kappa) {
    list(dMYZdt(t, y, pars, Lambda, kappa))
  }, parms = params, method = 'lsoda', Lambda = Lambda, kappa = kappa
  )

  # equilibrium solutions (forward)
  M_eq <- as.vector(Omega_inv %*% Lambda)
  M_sim <- as.vector(out[2, params$ix$MYZ$M_ix+1])

  P_eq <- as.vector(solve(diag(f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  P_sim <- as.vector(out[2, params$ix$MYZ$P_ix+1])

  Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Y_sim <- as.vector(out[2, params$ix$MYZ$Y_ix+1])

  Z_eq <- as.vector(Omega_inv %*% OmegaEIP %*% diag(f*q*kappa) %*% (M_eq - Y_eq))
  Z_sim <- as.vector(out[2, params$ix$MYZ$Z_ix+1])

  expect_equal(M_eq, M_sim, tolerance = numeric_tol)
  expect_equal(P_eq, P_sim, tolerance = numeric_tol)
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

library(expm)
library(deSolve)

test_that("RM models reach equilibrium", {
  nPatches <- 3
  f <- 0.3
  q <- 0.9
  g <- 1/20
  sigma <- 1/10
  tau <- 11
  nu <- 1/2
  eggsPerBatch <- 30

  K <- matrix(0, nPatches, nPatches)
  K[1, 2:3] <- c(0.2, 0.8)
  K[2, c(1,3)] <- c(0.5, 0.5)
  K[3, 1:2] <- c(0.7, 0.3)
  K <- t(K)

  Omega <- diag(g, nPatches) + ((diag(nPatches) - K) %*% diag(sigma, nPatches))
  OmegaEIP <- expm::expm(-Omega * tau)

  kappa <- c(0.1, 0.075, 0.025)
  Lambda <- c(5, 10, 8)

  params <- list(
    nPatches = nPatches,
    M_ix = 1:3,
    G_ix = 4:6,
    Y_ix = 7:9,
    Z_ix = 10:12
  )

  # ODE
  MYZpar <- make_parameters_MYZ_RM_ode(Omega = Omega, OmegaEIP = OmegaEIP, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = rep(0, nPatches), G0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 = rep(0, nPatches))
  params$MYZpar <- MYZpar

  y0 <- rep(0, 12)

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, Lambda, kappa) {
    list(dMYZdt(t, y, pars, Lambda, kappa))
  }, parms = params, method = 'lsoda', Lambda = Lambda, kappa = kappa)

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

  expect_true(all(approx_equal(M_eq, M_sim, tol = 1e-4)))
  expect_true(all(approx_equal(G_eq, G_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Y_eq, Y_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Z_eq, Z_sim, tol = 1e-4)))

  # equilibrium solutions (backward)
  MY_eq <- as.vector(diag(1/(f*q*kappa)) %*% OmegaEIP_inv %*% Omega %*% Z_eq)
  MY_sim <- M_sim - Y_sim

  Y_eq <- as.vector(Omega_inv %*% diag(f*q*kappa) %*% MY_eq)

  M_eq <- MY_eq + Y_eq

  Lambda_eq <- as.vector(Omega %*% M_eq)

  expect_true(all(approx_equal(MY_eq, MY_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Y_eq, Y_sim, tol = 1e-4)))
  expect_true(all(approx_equal(M_eq, M_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Lambda_eq, Lambda, tol = 1e-4)))

  # DDE
  MYZpar <- make_parameters_MYZ_RM_dde(Omega = Omega, OmegaEIP = OmegaEIP, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, tau = tau, M0 = rep(0, nPatches), G0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 = rep(0, nPatches))
  params$MYZpar <- MYZpar

  out <- deSolve::dede(y = y0, times = c(0, 365), func = function(t, y, pars, Lambda, kappa) {
    list(dMYZdt(t, y, pars, Lambda, kappa))
  }, parms = params, method = 'lsoda', Lambda = Lambda, kappa = rbind(kappa, kappa))

  # equilibrium solutions (forward)
  M_eq <- as.vector(Omega_inv %*% Lambda)
  M_sim <- as.vector(out[2, params$M_ix+1])

  G_eq <- as.vector(solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  G_sim <- as.vector(out[2, params$G_ix+1])

  Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Y_sim <- as.vector(out[2, params$Y_ix+1])

  Z_eq <- as.vector(Omega_inv %*% OmegaEIP %*% diag(f*q*kappa) %*% (M_eq - Y_eq))
  Z_sim <- as.vector(out[2, params$Z_ix+1])

  expect_true(all(approx_equal(M_eq, M_sim, tol = 1e-4)))
  expect_true(all(approx_equal(G_eq, G_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Y_eq, Y_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Z_eq, Z_sim, tol = 1e-4)))

  # equilibrium solutions (backward)
  MY_eq <- as.vector(diag(1/(f*q*kappa)) %*% OmegaEIP_inv %*% Omega %*% Z_eq)
  MY_sim <- M_sim - Y_sim

  Y_eq <- as.vector(Omega_inv %*% diag(f*q*kappa) %*% MY_eq)

  M_eq <- MY_eq + Y_eq

  Lambda_eq <- as.vector(Omega %*% M_eq)

  expect_true(all(approx_equal(MY_eq, MY_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Y_eq, Y_sim, tol = 1e-4)))
  expect_true(all(approx_equal(M_eq, M_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Lambda_eq, Lambda, tol = 1e-4)))
})

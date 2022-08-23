library(MASS)
library(expm)
library(deSolve)

test_that("forced emergence works with equilibrium", {
  nPatches <- 3
  nHabitats <- 4
  f <- 0.3
  q <- 0.9
  g <- 1/20
  sigma <- 1/10
  tau <- 11
  nu <- 1/2
  eggsPerBatch <- 30

  calN <- matrix(0, nPatches, nHabitats)
  calN[1,1] <- 1
  calN[2,2] <- 1
  calN[3,3] <- 1
  calN[3,4] <- 1

  calU <- matrix(0, nHabitats, nPatches)
  calU[1,1] <- 1
  calU[2,2] <- 1
  calU[3:4,3] <- 0.5

  K <- matrix(0, nPatches, nPatches)
  K[1, 2:3] <- c(0.2, 0.8)
  K[2, c(1,3)] <- c(0.5, 0.5)
  K[3, 1:2] <- c(0.7, 0.3)
  K <- t(K)

  Omega <- diag(g, nPatches) + ((diag(nPatches) - K) %*% diag(sigma, nPatches))
  OmegaEIP <- expm::expm(-Omega * tau)

  kappa <- c(0.1, 0.075, 0.025)
  Lambda <- c(5, 10, 8)

  # equilibrium solutions (forward)
  Omega_inv <- solve(Omega)
  OmegaEIP_inv <- expm::expm(Omega * tau)

  M_eq <- as.vector(Omega_inv %*% Lambda)
  G_eq <- as.vector(solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Z_eq <- as.vector(Omega_inv %*% OmegaEIP %*% diag(f*q*kappa) %*% (M_eq - Y_eq))

  # the "Lambda" for the dLdt model
  alpha <- as.vector(ginv(calN) %*% Lambda)

  params <- list(
    nPatches = nPatches,
    nHabitats = nHabitats,
    calU = calU,
    calN = calN,
    M_ix = 1:3,
    G_ix = 4:6,
    Y_ix = 7:9,
    Z_ix = 10:12
  )

  # ODE
  MYZpar <- make_parameters_MYZ_RM_ode(Omega = Omega, OmegaEIP = OmegaEIP, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = rep(0, nPatches), G0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 = rep(0, nPatches))
  Lpar <- make_parameters_L_trace(Lambda = alpha)
  params$MYZpar <- MYZpar
  params$Lpar <- Lpar

  y0 <- rep(0, 12)
  y0[params$M_ix] <- M_eq
  y0[params$G_ix] <- G_eq
  y0[params$Y_ix] <- Y_eq
  y0[params$Z_ix] <- Z_eq

  out <- deSolve::ode(y = y0, times = c(0, 365), func = xDE_diffeqn_mosy, parms = params, method = 'lsoda', kappa = kappa)

  M_sim <- as.vector(out[2, params$M_ix+1])
  G_sim <- as.vector(out[2, params$G_ix+1])
  Y_sim <- as.vector(out[2, params$Y_ix+1])
  Z_sim <- as.vector(out[2, params$Z_ix+1])

  expect_true(all(approx_equal(M_eq, M_sim, tol = 1e-4)))
  expect_true(all(approx_equal(G_eq, G_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Y_eq, Y_sim, tol = 1e-4)))
  expect_true(all(approx_equal(Z_eq, Z_sim, tol = 1e-4)))

})

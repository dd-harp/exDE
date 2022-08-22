library(deSolve)

test_that("human SIP model remains at equilibrium", {
  nStrata <- 3
  H <- c(100, 500, 250)
  X <- c(20, 120, 80)
  b <- 0.55
  c <- 0.15
  r <- 1/200
  eta <- c(1/30, 1/40, 1/35)
  rho <- c(0.05, 0.1, 0.15)

  P <- diag(1/eta) %*% diag(rho/(1-rho)) %*% (r*X)
  EIR <- diag(1/b, nStrata) %*% diag(1/(1-rho)) %*% ((r*X)/(H-X-P))

  params <- list(
    nStrata = nStrata,
    X_ix = 1:3,
    P_ix = 4:6
  )

  Xpar <- make_parameters_X_SIP(b = b, c = c, r = r, rho = rho, eta = eta, X0 = X, P0 = as.vector(P))
  params$Xpar <- Xpar

  y0 <- rep(0, 6)
  y0[params$X_ix] <- X
  y0[params$P_ix] <- P

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, EIR) {
    list(dXdt(t, y, pars, EIR))
  }, parms = params, method = 'lsoda', EIR = as.vector(EIR))

  expect_true(all(approx_equal(as.vector(out[2L, params$X_ix+1]), X)))
  expect_true(all(approx_equal(as.vector(out[2L, params$P_ix+1]), as.vector(P))))
})

library(deSolve)

numeric_tol <- 1e-5

test_that("human SIP model remains at equilibrium", {
  nStrata <- 3
  H <- c(100, 500, 250)
  membershipH = 1:nStrata
  searchWtsH = rep(1, nStrata)
  X <- c(20, 120, 80)
  b <- 0.55
  c <- 0.15
  r <- 1/200
  eta <- c(1/30, 1/40, 1/35)
  rho <- c(0.05, 0.1, 0.15)
  Psi <- matrix(data = 1,nrow = 1, ncol = nStrata)

  P <- diag(1/eta) %*% diag(rho/(1-rho)) %*% (r*X)
  EIR <- diag(1/b, nStrata) %*% diag(1/(1-rho)) %*% ((r*X)/(H-X-P))

  params <- make_parameters_xde()
  params$nStrata <- nStrata
  params$nPatches <- 1

  params = make_parameters_demography_null(pars = params, H=H, membershipH=membershipH,
                                           searchWtsH=searchWtsH, TimeSpent=Psi)
  params = make_parameters_X_SIP(pars = params, b = b, c = c, r = r, eta=eta, rho=rho)
  params = make_inits_X_SIP(pars = params, X, P)

  params = make_indices(params)


  # set initial conditions
  y0 <- get_inits(params)


  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, EIR) {
    list(dXdt(t, y, pars, EIR))
  }, parms = params, method = 'lsoda', EIR = as.vector(EIR))

  expect_equal(as.vector(out[2L, params$X_ix+1]), X, tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$P_ix+1]), as.vector(P), tolerance = numeric_tol)
})

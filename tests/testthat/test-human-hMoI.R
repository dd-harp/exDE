library(deSolve)

numeric_tol <- 1e-5

test_that("human hybrid MoI model remains at equilibrium", {
  nStrata <- 3
  H <- c(100, 500, 250)
  b <- 0.55
  c1 <- 0.05
  c2 <- 0.25
  r1 <- 1/250
  r2 <- 1/50
  Psi <- matrix(data = 1,nrow = 1, ncol = nStrata)

  m20 <- 1.5
  h <- r2*m20
  m10 <- h/r1

  EIR <- h/b

  params <- list(
    nStrata = nStrata
  )
  params <- list2env(params)

  make_parameters_X_hMoI(pars = params, b = b, c1 = c1, c2 = c2, r1 = r1, r2 = r2, Psi = Psi, m10 = m10, m20 = m20, H = H)
  make_indices(params)

  y0 <- rep(0, 6)
  y0[params$m1_ix] <- m10
  y0[params$m2_ix] <- m20

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, EIR) {
    list(dXdt(t, y, pars, EIR))
  }, parms = params, method = 'lsoda', EIR = as.vector(EIR))

  expect_equal(as.vector(out[2L, params$m1_ix+1]), rep(m10, nStrata), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$m2_ix+1]), rep(m20, nStrata), tolerance = numeric_tol)

})

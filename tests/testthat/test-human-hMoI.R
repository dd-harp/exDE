library(deSolve)

numeric_tol <- 1e-5

test_that("human hybrid MoI model remains at equilibrium", {
  nStrata <- 3
  H <- c(100, 500, 250)
  residence = 1:nStrata
  searchWtsH = rep(1, nStrata)
  b <- 0.55
  c1 <- 0.05
  c2 <- 0.25
  r1 <- 1/250
  r2 <- 1/50
  TaR <- matrix(data = 1,nrow = 1, ncol = nStrata)

  m20 <- 1.5
  h <- r2*m20
  m10 <- h/r1

  EIR <- h/b

  params <- make_parameters_xde()
  params$nStrata = nStrata
  params$nPatches = 1
  params = make_parameters_demography_null(pars = params, H=H, residence=residence,
                                           searchWts=searchWtsH, TaR=TaR)
  params = make_parameters_X_hMoI(pars = params, b = b, c1 = c1, c2 = c2, r1 = r1, r2 = r2)
  params = make_inits_X_hMoI(pars = params, m10 = rep(m10,nStrata), m20 = rep(m20,nStrata))

  params = make_indices(params)

  # set initial conditions
  y0 <- get_inits(params)

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, EIR) {
    list(dXdt(t, y, pars, EIR))
  }, parms = params, method = 'lsoda', EIR = as.vector(EIR))

  expect_equal(as.vector(out[2L, params$Xpar$m1_ix+1]), rep(m10, nStrata), tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$Xpar$m2_ix+1]), rep(m20, nStrata), tolerance = numeric_tol)

})

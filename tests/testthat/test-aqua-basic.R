library(deSolve)


test_that("basic competition stays at equilibrium", {

  numeric_tol <- 1e-5

  nHabitats <- 3
  alpha <- c(10, 50, 20)
  eggs_laid <- c(250, 500, 170)
  psi <- 1/10
  phi <- 1/12

  L <- alpha/psi
  theta <- (eggs_laid - psi*L - phi*L)/(L^2)

  params <- make_parameters_xde()
  params$nHabitats = nHabitats

  # ODE
  params = make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta=theta)
  params = make_inits_L_basic(pars = params, L0 = L)
  params = make_indices(params)
  params$eggs_laid = list()
  params$eggs_laid[[1]] = eggs_laid

  y0 <- rep(0, 3)

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, s) {
    list(dLdt(t, y, pars, s))
  }, parms = params, method = 'lsoda', s=1)

  expect_equal(as.vector(out[2L, 2:4]), L, tolerance = numeric_tol)
})

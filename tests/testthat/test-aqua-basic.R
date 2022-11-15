library(deSolve)

test_that("basic competition stays at equilibrium", {
  nHabitats <- 3
  alpha <- c(10, 50, 20)
  eta <- c(250, 500, 170)
  psi <- 1/10
  phi <- 1/12

  L <- alpha/psi
  theta <- (eta - psi*L - phi*L)/(L^2)

  params <- list(
    nHabitats = nHabitats
  )
  params <- list2env(params)

  # ODE
  make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta = theta, L0 = L)
  make_indices(params)

  y0 <- rep(0, 3)

  out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, eta) {
    list(dLdt(t, y, pars, eta))
  }, parms = params, method = 'lsoda', eta = eta)

  expect_true(all(approx_equal(as.vector(out[2L, 2:4]), L)))
})

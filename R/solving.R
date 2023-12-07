#' @title Solve a system of equations
#' @description This method dispatches on the type of `pars$xde`.
#' @param pars a [list] that defines a model
#' @param Tmax the last time point, run from 0...Tmax
#' @param dt the time interval for outputs
#' @return a [list]
#' @export
xde_solve = function(pars, Tmax=365, dt=1){
  UseMethod("xde_solve", pars$xde)
}

#' @title Solve for the steady state or stable orbit of a system of equations
#' @description This method dispatches on the type of `pars$xde`.
#' @param pars a [list] that defines a model
#' @param Ymax the number of years to burn-in
#' @return a [list]
#' @export
xde_stable_orbit = function(pars, Ymax=10){
  pars <- xde_solve(pars, Tmax = Ymax*365, dt=1)
  deout = tail(pars$orbits$deout, 365)
  deout[,1] = c(1:365)
  steady <- parse_deout(deout, pars)
  steady$terms = pars$orbits
  for(i in 1:length(steady$terms))
    steady$terms[[i]] = tail(steady$terms[[i]],365)
  pars$outputs$stable_orbits <- steady
  return(pars)
}

#' @title Solve for the steady state of a system of equations using [rootSolve::steady]
#' @description This method dispatches on the type of `pars$xde`
#' @param pars a [list] that defines a model
#' @return a [list]
#' @export
xde_steady = function(pars){
  y0 = get_inits(pars)
  pars1 <- dde2ode_MYZ(pars)
  rootSolve::steady(y=y0, func = xDE_diffeqn, parms = pars1)$y -> y_eq
  pars$outputs$steady = parse_deout_vec(y_eq, pars)
  return(pars)
}

#' @title Solve a system of equations as an ode
#' @description Implements [xde_solve] for ordinary differential equations
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.ode = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn, parms = pars, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}

#' @title Solve a system of equations as a dde
#' @description Implements [xde_solve] for delay differential equations
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.dde = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::dede(y = y0, times = tt, func = xDE_diffeqn, parms = pars, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}

#' @title Solve a system of equations for aquatic dynamics, forced by egg deposition, using xDE_diffeqn_aquatic
#' @description Implements [xde_solve] for aquatic dynaamic
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.aqua = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits_L(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn_aquatic, parms = pars, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}

#' @title Solve a system of delay equations for aquatic dynamics, forced by egg deposition, using xDE_diffeqn_aquatic
#' @description Implements [xde_solve] for aquatic dynaamic
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.aqua_dde = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits_L(pars)
  deSolve::dede(y = y0, times = tt, func = xDE_diffeqn_aquatic, parms = pars, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}

#' @title Solve a system of equations for mosquito ecology using xDE_diffeqn_mosy
#' @description Implements [xde_solve] for mosquito dynamics (no transmission)
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.mosy = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn_mosy, parms = pars, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}

#' @title Solve a system of delay differential equations for mosquito ecology using xDE_diffeqn_mosy
#' @description Implements [xde_solve] for mosquito dynamics (no transmission)
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.mosy_dde = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::dede(y = y0, times = tt, func = xDE_diffeqn_mosy, parms = pars, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}

#' @title Solve a system of equations with xDE_diffeqn_human
#' @description Implements [xde_solve] for mosquito dynamics (no transmission)
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.human = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn_human, parms = pars, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}


#' @title Solve a system of equations with xDE_diffeqn_cohort
#' @description Implements [xde_solve] for mosquito dynamics (no transmission)
#' @inheritParams xde_solve
#' @return a [list]
#' @export
xde_solve.cohort = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn_cohort, parms=pars, F_eir = pars$F_eir, method = "lsoda") -> out
  pars$outputs$orbits = parse_deout(out, pars)
  return(pars)
}

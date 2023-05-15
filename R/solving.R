#' @title Solve a system of equations
#' @description This method dispatches on the type of `pars$xde`.
#' @param pars a [list] that defines a model
#' @param Tmax the last time point, run from 0...Tmax
#' @param dt the time interval for outputs
#' @return a [matrix]
#' @export
xde_solve = function(pars, Tmax=365, dt=1){
  UseMethod("xde_solve", pars$xde)
}

#' @title Solve a system of equations as an ode
#' @description Implements [xde_solve] for ordinary differential equations
#' @inheritParams xde_solve
#' @return a [matrix]
#' @export
xde_solve.ode = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn, parms = pars, method = "lsoda")
}

#' @title Solve a system of equations as a dde
#' @description Implements [xde_solve] for delay differential equations
#' @inheritParams xde_solve
#' @return a [matrix]
#' @export
xde_solve.dde = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::dede(y = y0, times = tt, func = xDE_diffeqn, parms = pars, method = "lsoda")
}

#' @title Solve a system of equations for aquatic dynamics, forced by egg deposition, using xDE_diffeqn_aquatic
#' @description Implements [xde_solve] for aquatic dynaamic
#' @inheritParams xde_solve
#' @return a [matrix]
#' @export
xde_solve.aqua = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits_L(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn_aquatic, parms = pars, method = "lsoda")
}

#' @title Solve a system of equations for mosquito ecology using xDE_diffeqn_mosy
#' @description Implements [xde_solve] for mosquito dynamics (no transmission)
#' @inheritParams xde_solve
#' @return a [matrix]
#' @export
xde_solve.mosy = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn_mosy, parms = pars, method = "lsoda")
}

#' @title Solve a system of equations with xDE_diffeqn_human
#' @description Implements [xde_solve] for mosquito dynamics (no transmission)
#' @inheritParams xde_solve
#' @return a [matrix]
#' @export
xde_solve.human = function(pars, Tmax=365, dt=1){
  tt = seq(0, Tmax, by=dt)
  y0 = get_inits(pars)
  deSolve::ode(y = y0, times = tt, func = xDE_diffeqn_human, parms = pars, method = "lsoda")
}


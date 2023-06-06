
#' @title Compute the blood feeding rate, f
#' @description This method dispatches on the type of `pars$MYZpar$f_par`. It should
#' set the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f <- function(t, pars) {
  UseMethod("F_f", pars$MYZpar$f_par)
}

#' @title Static model for the blood feeding rate
#' @description Implements [F_f] for a static model
#' @inheritParams F_f
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f.static <- function(t, pars){
  pars$MYZpar$f0
}

#' @title Type 2 functional response for the blood feeding rate
#' @description Implements [F_f] for a static model
#' @inheritParams F_f
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f.type2 <- function(t, pars){
  B = pars$B
  with(pars$MYZpar$f_par,{
    return(fx*sf*B/(1+sf*B))
  })
}

#' @title Compute the human blood fraction
#' @description This method dispatches on the type of `pars$MYZpar$q_par`. It should
#' set the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q <- function(t, pars) {
  UseMethod("F_q", pars$MYZpar$q_par)
}

#' @title Static model for human blood fraction
#' @description Implements [F_q] for a static model
#' @inheritParams F_q
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q.static <- function(t, pars){
  pars$MYZpar$q0
}

#' @title Static model for human blood fraction
#' @description Implements [F_q] for a static model
#' @inheritParams F_q
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q.dynamic <- function(t, pars){
  with(pars,{
    return((W+Visitors)/B)
  })
}

#' @title Compute mosguito survival
#' @description This method dispatches on the type of `pars$MYZpar$g_par`. It should
#' set the values of g to (possibly changing) baseline values.
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_g <- function(t, pars) {
  UseMethod("F_g", pars$MYZpar$g_par)
}

#' @title Static model for mosquito survival
#' @description Implements [F_g] for a static model
#' @inheritParams F_g
#' @return a [numeric] vector of length `nPatches`
#' @export
F_g.static <- function(t, pars){
  pars$MYZpar$g0
}

#' @title Compute mosquito emigration rates
#' @description This method dispatches on the type of `pars$MYZpar$sigma_par`. It should
#' set the values of sigma to (possibly changing) baseline value(s).
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma <- function(t, pars) {
  UseMethod("F_sigma", pars$MYZpar$sigma_par)
}

#' @title Static model for mosquito emigration
#' @description Implements [F_sigma] for a static model
#' @inheritParams F_sigma
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma.static <- function(t, pars){
  pars$MYZpar$sigma0
}

#' @title Model for mosquito emigration based on resource availability
#' @description Implements [F_sigma] for a static model
#' @inheritParams F_sigma
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma.BQS <- function(t, pars){
  with(pars, with(MYZpar$sigma_par,{
    return(sigma_x*(sigma_B/(1+sB*B) + sigma_Q/(1+sQ*Q) + sigma_S/(1+sS*S)))
  }))
}


#' @title Compute the egg laying rate
#' @description This method dispatches on the type of `pars$MYZpar$nu_par`. It should
#' set the values of nu to (possibly changing) baseline value(s).
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu <- function(t, pars) {
  UseMethod("F_nu", pars$MYZpar$nu_par)
}

#' @title Static model for the egg laying rate
#' @description Implements [F_nu] for a static model
#' @inheritParams F_nu
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu.static <- function(t, pars){
  pars$MYZpar$nu0
}

#' @title Type 2 functional response for the blood feeding rate
#' @description Implements [F_nu] for a static model
#' @inheritParams F_nu
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu.type2 <- function(t, pars){
  Q = pars$Q
  with(pars$MYZpar$nu_par,{
    return(nux*snu*Q/(1+snu*Q))
  })
}

#' @title Compute the eip
#' @description This method dispatches on the type of `pars$MYZpar$eip_par`. It should
#' set the values of the eip
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eip <- function(t, pars) {
  UseMethod("F_eip", pars$MYZpar$eip_par)
}

#' @title Static model for human blood fraction
#' @description Implements [F_eip] for a static model
#' @inheritParams F_eip
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eip.static <- function(t, pars){
  pars$MYZpar$eip
}

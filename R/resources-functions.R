#' @title Set the availability of hosts
#' @description This method dispatches on the type of `pars$Wpar`. It should
#' compute the availability of the pathogen's hosts for blood feeding at time `t`
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
computeW <- function(t, y, pars) {
  UseMethod("computeW", pars$Wpar)
}

#' @title Static model for hosts
#' @description Implements [computeW] for a static model
#' @inheritParams computeW
#' @return a [numeric] vector of length `nPatches`
#' @export
computeW.static <- function(t, y, pars){
  pars$W
}

#' @title Static model for hosts
#' @description Implements [computeW] for a static model
#' @inheritParams computeW
#' @return a [numeric] vector of length `nPatches`
#' @export
computeW.dynamic<- function(t, y, pars){
  H = F_H(t, y, pars)
  pars$W = with(pars$Hpar, TaR %*% (wts_f*H))
}

#' @title Set the availability of hosts
#' @description This method dispatches on the type of `pars$Bpar`. It should
#' compute the availability of all blood hosts at time `t`
#' @param t current simulation time
#' @param pars an [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
computeB <- function(t, pars) {
  UseMethod("computeB", pars$Bpar)
}

#' @title Static model for hosts
#' @description Implements [computeB] for a static model
#' @inheritParams computeB
#' @return a [numeric] vector of length `nPatches`
#' @export
computeB.static <- function(t, pars){
  pars$B
}

#' @title Static model for hosts
#' @description Implements [computeB] for a dynamic model
#' @inheritParams computeB
#' @return a [numeric] vector of length `nPatches`
#' @export
computeB.dynamic <- function(t, pars){
  with(pars, return(W + Visitors + other^pars$Bpar$zeta))
}

#' @title Set the availability of aquatic habitats
#' @description This method dispatches on the type of `pars$Qpar`. It should
#' compute the availability of aquatic habitats at time `t`
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
computeQ <- function(t, pars) {
  UseMethod("computeQ", pars$Qpar)
}

#' @title Static model for hosts
#' @description Implements [computeQ] for a static model
#' @inheritParams computeQ
#' @return a [numeric] vector of length `nPatches`
#' @export
computeQ.static <- function(t, pars){
  pars$Q
}

#' @title Static model for hosts
#' @description Implements [computeQ] for a static model
#' @inheritParams computeQ
#' @return a [numeric] vector of length `nPatches`
#' @export
computeQ.dynamic <- function(t, pars){
  pars$Q = pars$calN %*% pars$Lpar$searchQ
}

#' @title Set the availability of other blood hosts
#' @description This method dispatches on the type of `pars$OBpar`. It should
#' compute the availability of other blood hosts at time `t`
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_other <- function(t, pars) {
  UseMethod("F_other", pars$OBpar)
}

#' @title Static model for other blood hosts
#' @description Implements [F_other] for a static model
#' @inheritParams F_other
#' @return a [numeric] vector of length `nPatches`
#' @export
F_other.static <- function(t, pars){
  pars$MYZpar$other
}

#' @title Set the availability of sugar
#' @description This method dispatches on the type of `pars$SGRpar`. It should
#' set the values of sugar availability at time `t`
#' @param t current simulation time
#' @param pars a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sugar <- function(t, pars) {
  UseMethod("F_sugar", pars$SGRpar)
}

#' @title Static model for sugar blood hosts
#' @description Implements [F_sugar] for a static model
#' @inheritParams F_sugar
#' @return a [numeric] vector of length `nStrata`
#' @export
F_sugar.static <- function(t, pars){
  pars$MYZpar$sugar
}

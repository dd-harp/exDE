# specialized methods for the aquatic mosquito trace model

#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the trace (forced emergence) model.
#' @inheritParams F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.trace <- function(t, y, pars) {
  with(pars$Lpar, Lambda*Lt(t, pars))
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the trace (forced emergence) model.
#' @inheritParams dLdt
#' @return a [numeric] vector
#' @export
dLdt.trace <- function(t, y, pars, eta) {
  numeric(0)
}

#' @title Setup Lpar.trace
#' @description Implements [setup_L] for the trace model
#' @inheritParams setup_L
#' @return a [list] vector
#' @export
setup_L.trace = function(pars, Lname,
                            membership=1, searchQ=1,
                            Lopts=list()){


  nHabitats = length(membership)
  pars$nHabits = nHabitats
  pars$calN = make_calN(pars$nPatches, membership)
  pars$calU = t(pars$calN)

  with(Lopts,{
    pars$Lname = "trace"
    pars = make_Lpar_trace(pars, Lopts)
    pars$Linits = numeric(0)
    return(pars)
})}

#' @title Make parameters for trace aquatic mosquito model
#' @param pars a [list]
#' @param Lopts a [list] that overwrites default values
#' @param Lambda vector of mean emergence rates from each aquatic habitat
#' @param Lt is a [function] of the form Lt(t,pars) that computes temporal fluctuations
#' @return none
#' @export
make_Lpar_trace = function(pars, Lopts=list(),
                           Lambda=1, Lt = NULL){
  with(Lopts,{
    Lpar = list()
    class(Lpar) <- "trace"
    Lpar$Lambda = checkIt(Lambda, pars$nHabitats)
    if(is.null(Lt)) Lt = function(t, pars){1}
    Lpar$Lt = Lt
    pars$Lpar <- Lpar
    return(pars)
})}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_L] for trace (forced emergence) model.
#' @inheritParams make_indices_L
#' @return none
#' @export
make_indices_L.trace <- function(pars) {
  pars$L_ix <- integer(0)
  return(pars)
}

#' @title Make parameters for trace aquatic mosquito model
#' @param pars a [list]
#' @param Lambda vector of mean emergence rates from each aquatic habitat
#' @param Lt is a [function] of the form Lt(t,pars) that computes temporal fluctuations
#' @return none
#' @export
make_parameters_L_trace <- function(pars, Lambda, Lt=function(t,pars){1}) {
  stopifnot(is.numeric(Lambda))
  Lpar <- list()
  class(Lpar) <- 'trace'
  Lpar$Lambda <- Lambda
  Lpar$Lt = Lt
  pars$Lpar <- Lpar
  return(pars)
}

#' @title Make parameters for trace aquatic mosquito model
#' @param pars a [list]
#' @param L0 is set to NULL for the trace model
#' @return none
#' @export
make_inits_L_trace<- function(pars, L0=NULL) {
  pars$Linits = numeric(0)
  return(pars)
}

#' @title Parse the variable names for the trace model
#' @description Implements [parse_deout_L] for the trace model.
#' @inheritParams parse_deout_L
#' @return varslist a [list]
#' @export
parse_deout_L.trace <- function(varslist, deout, pars) {
  return(varslist)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_L] for the GeRM model.
#' @inheritParams get_inits_L
#' @return none
#' @export
get_inits_L.trace <- function(pars){
  numeric(0)
}


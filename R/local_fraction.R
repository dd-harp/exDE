
#' @title Compute the local fraction
#' @description This method dispatches on the type of `pars$BFpar`
#' @param pars, a [list]
#' @return [list]
#' @export
LocalFrac <- function(pars){
  UseMethod("LocalFrac", pars$BFpar)
}

#' @title Compute transmission terms dynamically, no update required
#' @description This method dispatches on the type of `pars$BFpar`
#' @inheritParams LocalFrac
#' @return [list]
#' @export
LocalFrac.static <- function(pars){
  return(pars)
}

#' @title Compute transmission terms dynamically, no update required
#' @description This method dispatches on the type of `pars$BFpar`
#' @inheritParams LocalFrac
#' @return [list]
#' @export
LocalFrac.dynamic <- function(pars){
  pars = compute_local_frac(pars)
  return(pars)
}

#' @title Compute the local fraction
#' @description Compute the availability for the pathogen's hosts for blood feeding
#' @param pars a [list]
#' @return pars a [list]
#' @export
compute_local_frac <- function(pars){

  for(s in 1:pars$nVectors)
    pars$vars$local_frac[[s]] = with(pars$vars, W[[s]]/(W[[s]]+Visitors[[s]]))

  return(pars)
}

#' @title Set up the local_fraction for static models
#' @description The local fraction
#' @param pars a [list]
#' @param local_frac is the fraction of the ambient human / host population that is not a visitor
#' @return [list]
#' @export
setup_local_fraction_simple <- function(pars, local_frac=1) {

  pars$vars$local_frac = list()
  pars$vars$local_frac[[1]] = local_frac

  return(pars)
}


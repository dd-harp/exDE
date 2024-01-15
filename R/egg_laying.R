
#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
#' @export
EggLaying = function(t, y, pars){
  UseMethod("EggLaying", pars$EGGpar)
}

#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @inheritParams EggLaying
#' @return a [list]
#' @export
EggLaying.static = function(t, y, pars){
  pars = compute_EggsLaid(t, y, pars)
  return(pars)
}

#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @inheritParams EggLaying
#' @return a [list]
#' @export
EggLaying.simple = function(t, y, pars){
  for(s in 1:pars$nVectors)
    pars$eggs_laid[[s]] = pars$calU[[s]] %*% F_eggs(t, y, pars, s)
  return(pars)
}

#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @inheritParams EggLaying
#' @return a [list]
#' @export
EggLaying.setup = function(t, y, pars){
  for(s in 1:pars$nVectors) pars = compute_AvailableHabitat(pars, s)
  pars = compute_calU(t, y, pars)
  class(pars$EGGpar) <- "static"
  return(pars)
}


#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @inheritParams EggLaying
#' @return a [list]
#' @export
EggLaying.dynamic = function(t, y, pars){
  pars = compute_calU(t, y, pars)
  pars = compute_EggsLaid(t, y, pars)
  return(pars)
}

#' @title Compute calU
#' @description Computes the egg laying matrices
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
compute_calU = function(t, y, pars){

  for(s in 1:pars$nVectors)
    pars = make_calU_s(pars, s)
  return(pars)
}

#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
compute_EggsLaid = function(t, y, pars){
  for(s in 1:pars$nVectors)
    pars$eggs_laid[[s]] = pars$calU[[s]] %*% (F_eggs(t, y, pars, s)*pars$vars$Qfrac[[s]])

  return(pars)
}

#' @title Make the egg distribution matrix, calU
#' @param calN the habitat membership matrix
#' @param searchWts the habitat search weights
#' @return a [matrix] of dimensions `nHabitats` by `nPatches`
#' @export
make_calU = function(calN, searchWts=1){
  calU = searchWts*t(calN)
  colNorms = colSums(calU)
  ix = which(colNorms == 0)
  if(length(ix)>0) colNorms[ix]=1
  calU = calU %*% diag(1/colNorms)
  return(calU)
}

#' @title Make the egg distribution matrix, calU
#' @param pars the model object
#' @param s vector species index
#' @return [list]
#' @export
make_calU_s = function(pars, s){
  pars$calU[[s]] = make_calU(pars$calN, pars$EGGpar$searchWts[[s]])
  return(pars)
}


#' @title Setup the structures required for egg laying
#' @description This sets up the object EGGpar, which holds information
#' about the egg laying model. Several terms are expected to differ by
#' vector species, including the search weights (searchWts), the fraction
#' of eggs laid in habitat that is sometimes suitable (Qfrac), the
#' egg deposition matrix (calU), and the habitat egg laying rate (eggs_laid).
#' This function sets up lists.
#' @param pars the model object
#' @return a [list] vector
#' @export
setup_EGGpar_static = function(pars){
  up <- list()
  class(up) <- "static"

  pars$EGGpar <- up
  pars$EGGpar$searchWts = list()

  pars$vars$Qfrac = list()
  pars$vars$Qfrac[[1]] <- 1

  pars$calU = list()
  pars$eggs_laid = list()

  return(pars)
}

#' @title Setup egg laying for most models
#' @description Sets up the egg-deposition matrix calU for the s^th species
#' @param pars the model object
#' @param s the vector species index
#' @param searchQ the membership matrix
#' @param Lopts a [list] of options to override defaults
#' @return a [list] vector
#' @export
setup_EggLaying_static = function(pars, s, searchQ=1, Lopts=list()){with(Lopts,{

  # Habitat search weights
  searchQ = checkIt(searchQ, pars$nHabitats)
  pars$EGGpar$searchWts[[s]] = checkIt(searchQ, pars$nHabitats, F)
  pars <- make_calU_s(pars, s)
  return(pars)
})}

#' @title Setup egg laying for most models
#' @description Sets up the egg-deposition matrix calU for the s^th species
#' @param pars the model object
#' @param s the vector species index
#' @param searchQ the membership matrix
#' @param Lopts a [list] of options to override defaults
#' @return a [list] vector
#' @export
setup_EggLaying_simple = function(pars, s, searchQ=1, Lopts=list()){with(Lopts,{

  # Habitat search weights
  searchQ = checkIt(searchQ, pars$nHabitats)
  pars$EGGpar$searchWts[[s]] = searchQ
  pars <- make_calU_s(pars, 1)
  class(pars$EGGpar) <- "simple"
  return(pars)
})}

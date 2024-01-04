#' @title Setup egg laying for most models
#' @description Sets up the egg-deposition matrix calU for the s^th species
#' @param pars the model object
#' @param s the vector species index
#' @param searchQ the membership matrix
#' @param Lopts a [list] of options to override defaults
#' @return a [list] vector
#' @export
setup_EggLaying = function(pars, s, searchQ=1, Lopts=list()){with(Lopts,{
  pars$search4habitat[[s]] = checkIt(searchQ, pars$nHabitats, F)
  pars$calU[[s]] = make_calU(pars$calN, searchQ)
  return(pars)
})}

#' @title Compute eggs laid
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
#' @export
EggLaying = function(t, y, pars){
  for(s in 1:pars$nVectors)
    pars$eta[[s]] = pars$calU[[s]] %*% F_eggs(t, y, pars,s)
  return(pars)
}

#' @title Searching for Habitat
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
#' @export
HabitatSearch = function(t, y, pars){
  UseMethod("HabitatSearch", pars$calU)
}

#' @title Searching for Habitat
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
#' @export
HabitatSearch.static = function(t, y, pars){
  return(pars)
}

#' @title Searching for Habitat
#' @description Computes eggs laid for the s^th species
#' @param t the time
#' @param y the state variables
#' @param pars the model object
#' @return a [list]
#' @export
HabitatSearch.dynamic = function(t, y, pars){
  for(s in 1:pars$nVectors)
    pars$calU[[s]] = with(pars, make_calU(calN, egg_laying[[s]]$searchQ))
  return(pars)
}

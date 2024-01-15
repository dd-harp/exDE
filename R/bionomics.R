
#' @title Set bionomic parameter rates relative to baseline
#' @description This calls Mbionomics and Lbionmics for each species. This function
#' resets bionomic parameters to their pre-control baseline value, which can later be
#' modified by vector control. In some models, the pre-control baseline is computed in
#' here as a function of resource availability.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @return a [list]
#' @export
Bionomics <- function(t, y, pars){
  UseMethod("Bionomics", pars$BIONOMICS)
}

#' @title Set bionomic parameter rates relative to baseline
#' @description This calls Mbionomics and Lbionmics for each species
#' @inheritParams Bionomics
#' @return a [list]
#' @export
Bionomics.static <- function(t, y, pars){
  return(pars)
}

#' @title Set bionomic parameter rates relative to baseline
#' @description This calls Mbionomics and Lbionmics for each species
#' @inheritParams Bionomics
#' @return a [list]
#' @export
Bionomics.setup <- function(t, y, pars){
 for(s in 1:pars$nVectors){
    pars <- MBionomics(t, y, pars, s)
    pars <- LBionomics(t, y, pars, s)
  }
  class(pars$BIONOMICS) <- "static"
  return(pars)
}

#' @title Set bionomic parameter rates relative to baseline
#' @description This calls Mbionomics and Lbionmics for each species
#' @inheritParams Bionomics
#' @return a [list]
#' @export
Bionomics.dynamic <- function(t, y, pars){
 for(s in 1:pars$nVectors){
    pars <- MBionomics(t, y, pars, s)
    pars <- LBionomics(t, y, pars, s)
  }
  return(pars)
}

#' @title Set bionomic parameter rates relative to baseline
#' @description This calls Mbionomics and Lbionmics for each species
#' @inheritParams Bionomics
#' @return a [list]
#' @export
Bionomics.Mdynamic <- function(t, y, pars){
  for(s in 1:pars$nVectors) pars <- MBionomics(t, y, pars, s)
  return(pars)
}

#' @title Set bionomic parameter rates relative to baseline
#' @description This calls Lbionomics for each species each time step
#' @inheritParams Bionomics
#' @return a [list]
#' @export
Bionomics.Ldynamic <- function(t, y, pars){
  for(s in 1:pars$nVectors) pars <- LBionomics(t, y, pars, s)
  return(pars)
}

#' @title Make parameters for the static model bionomics
#' @param pars a [list]
#' @return [list]
#' @export
setup_bionomics_static <- function(pars) {

  BIONOMICS <- list()
  class(BIONOMICS) <- "setup"
  pars$BIONOMICS <- BIONOMICS

  return(pars)
}

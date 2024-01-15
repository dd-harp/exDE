
#' @title Set up the blood feeding model
#' @description This sets up a list that stores all the information required by the
#' blood feeding model, including information about humans (the residence vector, search
#' weights for blood feeding mosquitoes); the mosquitoes circadian pattern (F_circadian);
#' and derived structures, including the time spent (TimeSpent) and the time at risk (TaR)
#' matrices. This also sets up the objects required for blood feeding, including host availability
#' for each species (Wi) and for all hosts (W), other blood hosts (Other), and total blood
#' availability (B). Finally, it sets up the lists that will hold the mixing matrix (beta), and
#' the transmission terms describing the entomological inoculation rate (EIR and eir), and the force
#' of infection (FoI), and the net infectiousness (kappa and ni).
#' @param pars a [list]
#' @return none
#' @export
setup_BFpar_static <- function(pars){
  up = list()
  class(up) <- "setup"

  # Human
  up$residence = list()
  up$searchWts = list()
  up$searchWts[[1]] = list()

  # Mosquito
  up$F_circadian = list()

  # Time Spent / Time at Risk
  up$TimeSpent = list()
  up$TaR = list()
  up$TaR[[1]] = list()

  # Mosquito
  pars$BFpar <- up

  # Available Blood Hosts, by vector species
  pars$vars$Wi = list()
  pars$vars$Wi[[1]] = list()
  pars$vars$W = list()
  pars$vars$Other = list()
  pars$vars$B = list()

  # Transmission Terms
  pars$beta = list()
  pars$beta[[1]] = list()

  pars$eir = list()
  pars$eir[[1]] = list()
  pars$EIR = list()
  pars$FoI = list()

  pars$ni = list()
  pars$ni[[1]] = list()
  pars$kappa = list()

  return(pars)
}

#' @title Set up blood feeding
#' @description This sets up a list that stores all the information
#' @param pars a [list]
#' @param i the host species index
#' @param s the vector species index
#' @param BFopts a [list]
#' @param residence is the patch where each stratum resides
#' @param searchWts is the blood feeding search weight for each stratum
#' @param F_circadian is a function that computes relative mosquito blood feeding activity rates by time of day
#' @return none
#' @export
setup_BloodFeeding <- function(pars, i, s=1, BFopts = list(), residence=1, searchWts=1, F_circadian=NULL){
  pars$BFpar$searchWts[[i]][[s]] = checkIt(searchWts, pars$Hpar[[i]]$nStrata)
  pars$BFpar$relativeBitingRate[[i]][[s]] = checkIt(searchWts, pars$Hpar[[i]]$nStrata)
  pars$BFpar$residence[[i]] = checkIt(residence, pars$Hpar[[i]]$nStrata)
  if(is.null(F_circadian)) pars$BFpar$F_circadian[[s]] = function(d){return(1)}
  return(pars)
}

#' @title Make TaR
#' @description Make a time at risk matrix (TaR) from a time spent matrix and a circadian function
#' @param t the time
#' @param pars a [list]
#' @param i the host species index
#' @param s the vector species index
#' @return none
#' @export
make_TaR <- function(t, pars, i, s){
  pars$BFpar$TaR[[i]][[s]] = with(pars$BFpar, TimeSpent[[i]]*F_circadian[[s]](t))
  return(pars)
}



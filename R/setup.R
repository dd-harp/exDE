# functions to set up models

#' @title Set up a model for xde_diffeqn
#' @param modelName is a name for the model (arbitrary)
#' @param MYZname is a character string defining a MYZ model
#' @param Xname is a character string defining a X model
#' @param Lname is a character string defining a L model
#' @param nPatches is the number of patches
#' @param HPop is the number of humans in each patch
#' @param membership is a vector that describes the patch where each aquatic habitat is found
#' @param MYZopts a list to configure the MYZ model
#' @param calK is either a calK matrix or a string that defines how to set it up
#' @param calKopts are the options to setup calK
#' @param Xopts a list to configure the X model
#' @param Hopts a list to configure the H model
#' @param residence is a vector that describes the patch where each human stratum lives
#' @param searchB is a vector of search weights for blood feeding
#' @param TaR is either a TaR matrix or a string to call a function that sets it up
#' @param TaRopts are the options to setup TaR
#' @param searchQ is a vector of search weights for egg laying
#' @param Lopts a list to configure the L model
#' @return a [list]
#' @export
xde_setup = function(modelName,

                     # Dynamical Components
                     MYZname = "RM",
                     Xname = "SIS",
                     Lname = "trace",

                     # Model Structure
                     nPatches = 1,
                     HPop=1000,
                     membership=1,

                     # Adult Mosquito Options
                     MYZopts = list(),
                     calK ="herethere",
                     calKopts = list(),

                     # Human Strata / Options
                     Xopts = list(),
                     Hopts = list(),
                     residence=1,
                     searchB = 1,
                     TaR = "athome",
                     TaRopts=list(),

                     # Aquatic Mosquito Options
                     searchQ = 1,
                     Lopts = list()

){

  pars = make_parameters_xde()
  pars$modelName = modelName

  # Structure
  pars$nPatches = nPatches
  pars$nStrata = length(HPop)
  pars$nHabitats = length(membership)

  pars = setup_Hpar(pars, HPop, residence, searchB, Hopts)
  pars$Hpar$TaR = make_TaR(nPatches, pars$Hpar$residence, TaR, TaRopts)

  # Dynamics
  calK = make_calK(nPatches, calK, calKopts)
  pars = setup_MYZ(pars, MYZname, nPatches, MYZopts, calK)
  pars = setup_L(pars, Lname, membership, searchQ, Lopts)
  pars = setup_X(pars, Xname, Xopts)

  pars = make_indices(pars)

  return(pars)
}

#' @title Set up a model for xde_diffeqn_mosy
#' @param modelName is a name for the model (arbitrary)
#' @param MYZname is a character string defining a MYZ model
#' @param Lname is a character string defining a L model
#' @param nPatches is the number of patches
#' @param membership is a vector that describes the patch where each aquatic habitat is found
#' @param MYZopts a list to configure the MYZ model
#' @param calK is either a calK matrix or a string that defines how to set it up
#' @param calKopts are the options to setup calK
#' @param searchQ is a vector of search weights for egg laying
#' @param Lopts a list to configure the L model
#' @param kappa values -- net infectivity to force adult infection dynamics
#' @return a [list]
#' @export
xde_setup_mosy = function(modelName,

                     # Dynamical Components
                     MYZname = "basicM",
                     Lname = "basic",

                     # Model Structure
                     nPatches = 1,
                     membership=1,

                     # Adult Mosquito Options
                     MYZopts = list(),
                     calK ="herethere",
                     calKopts = list(),

                     # Aquatic Mosquito Options
                     searchQ = 1,
                     Lopts = list(),

                     # forcing
                     kappa=NULL

){

  pars = make_parameters_xde()
  class(pars$xde) <- "mosy"
  pars$modelName = modelName

  # Structure
  pars$nPatches = nPatches
  pars$nHabitats = length(membership)

  # Dynamics
  calK = make_calK(nPatches, calK, calKopts)
  pars = setup_MYZ(pars, MYZname, nPatches, MYZopts, calK)
  pars = setup_L(pars, Lname, membership, searchQ, Lopts)

  if(is.null(kappa))  kappa = rep(0, nPatches)
  pars$kappa = checkIt(kappa, nPatches)

  pars$max_ix <- 0
  pars = make_indices_L(pars)
  pars = make_indices_MYZ(pars)

  return(pars)
}


#' @title Set up a model for xde_diffeqn_aqua
#' @param modelName is a name for the model (arbitrary)
#' @param nHabitats is the number of habitats
#' @param Lname is a character string defining a L model
#' @param Lopts a list to configure the L model
#' @param MYZopts a list to configure F_eggs from the Gtrace model
#' @param LSMname is a character string defining a LSM model
#' @return a [list]
#' @export
xde_setup_aquatic = function(modelName,
                     nHabitats = 1,
                     Lname = "basic",
                     Lopts = list(),
                     MYZopts = list(),
                     LSMname = "null"){

  pars = make_parameters_xde()
  class(pars$xde) <- "aqua"
  pars$modelName = modelName
  pars = setup_MYZ(pars, "Gtrace", nHabitats, MYZopts, NULL)
  pars$nHabitats = nHabitats
  membership = 1:nHabitats
  searchQ = rep(1, nHabitats)
  pars = setup_L(pars, Lname, membership, searchQ, Lopts)
  pars <- setup_lsm_null(pars)

  pars = make_indices(pars)

  return(pars)
}


#' @title Set up a model for xde_diffeqn_human
#' @param modelName is a name for the model (arbitrary)
#' @param Xname is a character string defining a X model
#' @param HPop is the number of humans in each patch
#' @param MYZopts a list to configure the MYZ model
#' @param Xopts a list to configure the X model
#' @param Hopts a list to configure the H model
#' @param residence is a vector that describes the patch where each human stratum lives
#' @param searchB is a vector of search weights for blood feeding
#' @param TaR is either a TaR matrix or a string to call a function that sets it up
#' @param TaRopts are the options to setup TaR
#' @return a [list]
#' @export
xde_setup_human = function(modelName,

                     # Dynamical Components
                     Xname = "SIS",

                     # Model Structure
                     HPop=1000,

                     # Adult Mosquito Options
                     MYZopts = list(),

                     # Human Strata / Options
                     Xopts = list(),
                     Hopts = list(),
                     residence=1,
                     searchB = 1,
                     TaR = "athome",
                     TaRopts=list()

){

  pars = make_parameters_xde()
  class(pars$xde) <- "human"
  pars$modelName = modelName

  # Structure
  nStrata = length(HPop)
  pars$nPatches = as.integer(nStrata)
  pars$nStrata = nStrata

  pars = setup_Hpar(pars, HPop, 1:nStrata, rep(1, nStrata), Hopts)
  pars$Hpar$TaR = make_TaR(pars$nPatches, pars$Hpar$residence, TaR, TaRopts)

  # Dynamics
  pars = setup_MYZ(pars, "Ztrace", pars$nPatches, MYZopts, calK=NULL)
  pars = setup_X(pars, Xname, Xopts)

  pars = make_indices(pars)

  return(pars)
}

#' @title Set up a model for xde_diffeqn_cohort
#' @param modelName is a name for the model (arbitrary)
#' @param F_eir is a function F_eir(t, pars) that returns the daily FoI
#' @param Xname is a character string defining a X model
#' @param HPop is the number of humans in each patch
#' @param Xopts a list to configure the X model
#' @param Hopts a list to configure the H model
#' @return a [list]
#' @export
xde_setup_cohort = function(modelName, F_eir,

                           # Dynamical Components
                           Xname = "SIS",

                           # Model Structure
                           HPop=1000,

                           # Human Strata / Options
                           Xopts = list(),
                           Hopts = list()

){

  pars = make_parameters_xde()
  class(pars$xde) <- "cohort"
  pars$modelName = modelName
  pars$F_eir = F_eir

  # Structure
  nStrata = length(HPop)
  pars$nPatches = as.integer(nStrata)
  pars$nStrata = nStrata

  pars = setup_Hpar(pars, HPop, 1:nStrata, rep(1, nStrata), Hopts)

  # Dynamics
  pars = setup_X(pars, Xname, Xopts)

  pars = make_indices(pars)

  return(pars)
}

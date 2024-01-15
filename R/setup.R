# functions to set up models

#' @title Set up a model for xde_diffeqn
#' @param modelName is a name for the model (arbitrary)
#' @param MYZname is a character string defining a MYZ model
#' @param Xname is a character string defining a X model
#' @param Lname is a character string defining a L model
#' @param nPatches is the number of patches
#' @param nVectors is the number of vector species
#' @param nHosts is the number of vertebrate host species
#' @param HPop is the number of humans in each patch
#' @param membership is a vector that describes the patch where each aquatic habitat is found
#' @param MYZopts a list to configure the MYZ model
#' @param calK is either a calK matrix or a string that defines how to set it up
#' @param calKopts are the options to setup calK
#' @param EIPname are the options to setup EIPmod
#' @param EIPopts are the options to setup EIPmod
#' @param Xopts a list to configure the X model
#' @param BFopts a list to configure the blood feeding model
#' @param residence is a vector that describes the patch where each human stratum lives
#' @param searchB is a vector of search weights for blood feeding
#' @param F_circadian is a function describing mosquito daily activity
#' @param TimeSpent is either a TimeSpent matrix or a string to call a function that sets it up
#' @param TimeSpentOpts are the options to setup TimeSpent
#' @param searchQ is a vector of search weights for egg laying
#' @param Lopts a list to configure the L model
#' @return a [list]
#' @export
xde_setup = function(modelName = "unnamed",

                     # Dynamical Components
                     MYZname = "RM",
                     Xname = "SIS",
                     Lname = "trace",

                     # Model Structure
                     nPatches = 1,
                     nVectors = 1,
                     nHosts = 1,
                     HPop=1000,
                     membership=1,

                     # Adult Mosquito Options
                     MYZopts = list(),
                     calK ="herethere",
                     calKopts = list(),
                     EIPname = "static",
                     EIPopts = list(),

                     # Human Strata / Options
                     Xopts = list(),

                     # Blood Feeding
                     BFopts = list(),
                     residence=1,
                     searchB = 1,
                     F_circadian = NULL,
                     TimeSpent = "athome",
                     TimeSpentOpts=list(),

                     # Aquatic Mosquito Options
                     searchQ = 1,
                     Lopts = list()

){

  pars = make_parameters_xde()
  class(pars$compute) = "xde"

  pars$modelName = modelName
  pars$Xname = Xname
  pars$MYZname = MYZname
  pars$Lname = Lname

  # Fixed Structural Elements
  pars$nPatches = nPatches
  pars$nVectors = nVectors
  pars$nHosts = nHosts
  pars$nHabitats = length(membership)
  pars$membership = membership
  pars$calN = make_calN(pars$nPatches, pars$membership)

  # Human Demography
  pars = setup_Hpar_static(pars, 1, HPop)

  # Blood Feeding
  pars = setup_BloodFeeding(pars, 1, 1, BFopts, residence, searchB, F_circadian)
  pars = make_TimeSpent(pars, 1, TimeSpent, TimeSpentOpts)

  # Adult Mosquito Dynamics
  EIPmod = setup_EIP(EIPname, EIPopts)
  calK = make_calK(nPatches, calK, calKopts)
  pars = setup_MYZpar(MYZname, pars, 1, MYZopts, EIPmod, calK)
  pars = setup_MYZinits(pars, 1, MYZopts)

  # Aquatic Mosquito Dynamics
  pars = setup_Lpar(Lname, pars, 1, Lopts)
  pars = setup_Linits(pars, 1, Lopts)

  pars = setup_EggLaying_static(pars, 1, searchQ)

  # Vertebrate Host Dynamics
  pars = setup_Xpar(Xname, pars, 1, Xopts)
  pars = setup_Xinits(pars, 1, Xopts)

  pars = make_indices(pars)

  y0 <- get_inits(pars)
  pars <- EggLaying(0, y0, pars)
  pars <- Resources(0, y0, pars)
  pars <- Bionomics(0, y0, pars)
  pars <- Transmission(0, y0, pars)

  return(pars)
}

#' @title Set up a model for xde_diffeqn_mosy
#' @param modelName is a name for the model (arbitrary)
#' @param MYZname is a character string defining a MYZ model
#' @param Lname is a character string defining a L model
#' @param nPatches is the number of patches
#' @param nVectors is the number of vector species
#' @param membership is a vector that describes the patch where each aquatic habitat is found
#' @param MYZopts a list to configure the MYZ model
#' @param calK is either a calK matrix or a string that defines how to set it up
#' @param calKopts are the options to setup calK
#' @param searchQ is a vector of search weights for egg laying
#' @param Lopts a list to configure the L model
#' @param kappa values -- net infectivity to force adult infection dynamics
#' @return a [list]
#' @export
xde_setup_mosy = function(modelName = "unnamed",

                     # Dynamical Components
                     MYZname = "basicM",
                     Lname = "basic",

                     # Model Structure
                     nPatches = 1,
                     nVectors = 1,
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
  class(pars$compute) = "na"

  pars$modelName = modelName
  pars$MYZname = MYZname
  pars$Lname = Lname

  # Structure
  pars$nPatches = nPatches
  pars$nHabitats = length(membership)
  pars$membership = membership
  pars$calN = make_calN(pars$nPatches, pars$membership)
  pars$nVectors = 1

  # Dynamics
  calK = make_calK(nPatches, calK, calKopts)
  pars = setup_MYZpar(MYZname, pars, 1, MYZopts, NULL, calK)
  pars = setup_MYZinits(pars, 1, MYZopts)

  # Aquatic Mosquito Dynamics
  pars = setup_Lpar(Lname, pars, 1, Lopts)
  pars = setup_Linits(pars, 1, Lopts)
  pars = setup_EggLaying_simple(pars, 1, searchQ)

  if(is.null(kappa))  kappa = rep(0, nPatches)
  pars$kappa[[1]] = checkIt(kappa, nPatches)

  pars = make_indices(pars)

  return(pars)
}


#' @title Set up a model for xde_diffeqn_aqua
#' @param modelName is a name for the model (arbitrary)
#' @param nHabitats is the number of habitats
#' @param nVectors is the number of vector species
#' @param Lname is a character string defining a L model
#' @param Lopts a list to configure the L model
#' @param MYZopts a list to configure F_eggs from the Gtrace model
#' @param LSMname is a character string defining a LSM model
#' @return a [list]
#' @export
xde_setup_aquatic = function(modelName = "unnamed",
                     nHabitats = 1,
                     nVectors = 1,
                     Lname = "basic",
                     Lopts = list(),
                     MYZopts = list(),
                     LSMname = "null"){

  pars = make_parameters_xde()
  class(pars$xde) <- "aqua"
  class(pars$compute) = "na"

  pars$modelName = modelName
  pars$MYZname = "Gtrace"
  pars$Lname = Lname

  pars$nVectors = nVectors
  pars = setup_MYZpar("Gtrace", pars, 1, MYZopts, EIPmod=NULL, calK=NULL)

  pars$nHabitats = nHabitats
  membership = 1:nHabitats
  pars$membership = membership
  pars$calN = make_calN(pars$nHabitats, pars$membership)
  searchQ = rep(1, nHabitats)
  pars = setup_Lpar(Lname, pars, 1, Lopts)
  pars = setup_Linits(pars, 1, Lopts)

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
#' @param BFopts a list to configure the blood feeding model
#' @param residence is a vector that describes the patch where each human stratum lives
#' @param searchB is a vector of search weights for blood feeding
#' @param F_circadian is a function describing mosquito daily activity
#' @param TimeSpent is either a TimeSpent matrix or a string to call a function that sets it up
#' @param TimeSpentOpts are the options to setup TimeSpent
#' @return a [list]
#' @export
xde_setup_human = function(modelName = "unnamed",

                     # Dynamical Components
                     Xname = "SIS",

                     # Model Structure
                     HPop=1000,

                     # Adult Mosquito Options
                     MYZopts = list(),

                     # Human Strata / Options
                     Xopts = list(),

                     # Blood Feeding
                     BFopts = list(),
                     residence=1,
                     searchB = 1,
                     F_circadian = NULL,
                     TimeSpent = "athome",
                     TimeSpentOpts=list()

){

  pars = make_parameters_xde()
  class(pars$xde) <- "human"
  class(pars$compute) = "human"

  pars$modelName = modelName
  pars$Xname = Xname
  pars$MYZname = "Ztrace"

  # Structure
  nStrata = length(HPop)
  pars$nPatches = as.integer(nStrata)
  pars$nStrata = nStrata

  pars = setup_Hpar_static(pars, 1, HPop)
  pars = setup_BloodFeeding(pars, 1, 1, BFopts, residence, searchB, F_circadian)
  pars = make_TimeSpent(pars, 1, TimeSpent, TimeSpentOpts)

  # Dynamics
  pars = setup_MYZpar("Ztrace", pars, 1, MYZopts, EIPmod=NULL, calK=NULL)

  pars = setup_Xpar(Xname, pars, 1, Xopts)
  pars = setup_Xinits(pars, 1, Xopts)

  pars = make_indices(pars)



  return(pars)
}

#' @title Set up a model for xde_diffeqn_cohort
#' @param F_eir is a function F_eir(t, pars) that returns the daily FoI
#' @param modelName is a name for the model (arbitrary)
#' @param Xname is a character string defining a X model
#' @param HPop is the number of humans in each patch
#' @param searchB is a vector of search weights for blood feeding
#' @param residence is a vector that describes the patch where each human stratum lives
#' @param Xopts a list to configure the X model
#' @return a [list]
#' @export
xde_setup_cohort = function(F_eir,
                           modelName = "unnamed",

                           # Dynamical Components
                           Xname = "SIS",

                           # Model Structure
                           HPop=1000,
                           searchB = 1,
                           residence = 1,

                           # Human Strata / Options
                           Xopts = list()

){

  pars = make_parameters_xde()
  class(pars$xde) <- "cohort"
  class(pars$compute) = "cohort"

  pars$modelName = modelName
  pars$Xname = Xname

  pars$F_eir = F_eir

  # Structure
  nStrata = length(HPop)
  pars$nPatches = as.integer(nStrata)
  pars$nHosts = 1

  pars = setup_Hpar_static(pars, 1, HPop)
  pars = setup_BloodFeeding(pars, 1, 1, list(), residence, searchB, NULL)

  # Dynamics
  pars = setup_Xpar(Xname, pars, 1, Xopts)
  pars = setup_Xinits(pars, 1, Xopts)

  pars = make_indices(pars)

  return(pars)
}

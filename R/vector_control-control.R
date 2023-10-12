
# specialized methods for the control model of vector control

#' @title Distribute vector control, the null model
#' @description Implements [VectorControl] for the control model of vector control (do nothing)
#' @inheritParams VectorControl
#' @return a named [list]
#' @export
VectorControl.control <- function(t, y, pars) {
  pars = DistributeBedNets(t, pars)
  pars = OwnBedNet(t, y, pars)
  pars = AreaSpray(t, pars)
  pars = SprayHouses(t, pars)
  pars = SugarBaits(t, pars)
  pars = TreatHabitats(t, pars)
  return(pars)
}

#' @title Vector control durability and effects
#' @description Implements [VectorControlEffects] for the control model of vector control (do nothing)
#' @inheritParams VectorControlEffects
#' @return a named [list]
#' @export
VectorControlEffects.control <- function(t, y, pars) {
  pars = BedNetEffects(t, pars)
  pars = AreaSprayEffects(t, pars)
  pars = IRSeffects(t, pars)
  pars = SugarBaitEffects(t, pars)
  pars = LSMeffects(t, pars)
  return(pars)
}

#' @title Distribute vector control, the null model
#' @description Implements [VectorControlEffectSizes] for the control model of vector control (do nothing)
#' @inheritParams VectorControlEffectSizes
#' @return a named [list]
#' @export
VectorControlEffectSizes.control <- function(t, y, pars) {
  pars = BedNetEffectSizes(t, pars)
  pars = AreaSprayEffectSizes(t, pars)
  pars = IRSeffectSizes(t, pars)
  pars = SugarBaitEffectSizes(t, pars)
  pars = LSMeffectSizes(t, pars)
  return(pars)
}

#' @title Make parameters for the control model of vector control (do nothing)
#' @param pars a [list]
#' @return none
#' @export
setup_vc_control <- function(pars) {
  pars = setup_control(pars)
  class(pars$VECTOR_CONTROL) <- 'control'
  pars <- setup_itn_null(pars)
  pars <- setup_area_spray_null(pars)
  pars <- setup_irs_null(pars)
  pars <- setup_sugar_baits_null(pars)
  pars <- setup_lsm_null(pars)
  return(pars)
}

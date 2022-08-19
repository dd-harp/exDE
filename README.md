# xDE (Extensible Differential Equations for mosquito-borne pathogen modeling)

<!-- badges: start -->
[![R-CMD-check](https://github.com/dd-harp/xDE/workflows/R-CMD-check/badge.svg)](https://github.com/dd-harp/xDE/actions)
<!-- badges: end -->
  
xDE provides tools to set up modular ordinary and delay differential equation spatial 
models for mosquito-borne pathogens, focusing on malaria. Modularity is achieved
by S3 dispatch on parameter lists for each component which is used to compute
the full set of differential equations. It can be regarded as the continuous-time
companion to the discrete stochastic [Micro-MoB](https://github.com/dd-harp/MicroMoB/tree/main)
framework.

## Modular Dynamics

Models for malaria transmission dynamics are naturally modular, structured by
vector life stage, host population strata, and by the spatial locations (patches) at
which transmission occurs (see figure below).

<p align="center">
  <img src="man/figures/modularity.png"/>
</p>

Two components describe mosquito ecology: dynamics of immature mosquitoes in aquatic habitats $dL/dt$ (blue component); and dynamics of adult mosquitoes $dM/dt$ (green components). Adults lay eggs ($\nu$) that are allocated among the habitats ($\eta$). The matrix $\mathcal{N}$ describes the spatial arrangement of patches and aquatic habitats. Some eggs will develop and emerge as mature adults $\alpha$, which are added to the mosquito population at each patch $\Lambda$.

Two additional components describe parasite infection and transmission (red): parasite infection dynamics in mosquitoes, $dY/dt$ and parasite infection dynamics in humans (purple), described by $dX/dt$. The biting distribution matrix $\beta$ links biting from mosquitoes in each patch to the available human population from each strata, which depends on how people spend their time across the landscape. The density of infectious mosquitoes $Z$ is used to compute the force of infection on humans, $h$; likewise the infectious human population $X$ is used to compute the net infectiousness $\kappa$ of humans on mosquitoes.

The interactions among these modules take place for a stratified human population within a spatial domain structured into patches that contain the aquatic habitats.

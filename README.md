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

### Adult Mosquitos

The generic functions which must be implemented for any specific model of adult mosquito dynamics are:

  * `F_EIR`:  compute $EIR$, which is a component in the force of infection $h$ upon humans.
  * `F_kappa`: computes $\kappa$, the net infectiousness of the human population to mosquitos.
  * `F_eggs`: compute $\nu$, the amount of eggs laid by mosquitos in each patch.
  * `dMYZdt`: compute the derivatives for all adult mosquitos (both uninfected $M$ and infected $Y$; $Z$ refers for infected _and_ infectious mosquitos, for models which make that distinction).

### Aquatic Mosquitos

The generic functions which must be implemented for any specific model of aquatic mosquito dynamics are:

  * `F_alpha`: compute $\alpha$, the rate of emergence of new adults from each aquatic habitat.
  * `dLdt`: compute the derivatives for aquatic mosquitos.
  
### Human Population

The generic functions which must be implemented for any specific model of human population dynamics are:

  * `F_x`: compute $X$, the net infectious population size.
  * `F_x_tau`: compute $X_{\tau}$, the lagged net infectious population size.
  * `dXdt`: compute the derivatives of the human population.

## Generalized equations

Because the framework does not make any assumptions on the specific internal dynamics
of each of the main components, only that they are able to provide the needed quantities
for the other components, a generalized set of differential equations can be written
to describe the dynamics of the system. Specific models can be used for each component
as needed. In `xDE` we use R's S3 dispatch to write generic differential equations
which can be specialized to specific models as needed.

The function `xDE_diffeqn` implements this generic differential equation model, which
has the following mathematical structure. This function is also responsible for computing
$\Lambda = \mathcal{U}\cdot\alpha$ and $\eta = \mathcal{N}\cdot\alpha$.
This structure is closely followed in the code.

$$
d{\mathcal{L}}/dt = F_{\mathcal{L}} \left(\eta, {\mathcal{L}} \right)
$$

$$
d {\mathcal{M}}/dt = F_{\mathcal{M}} \left(\Lambda, {\mathcal{M}} \right)
$$

$$
d {\mathcal{Y}}/dt = F_{\mathcal{Y}} \left(\kappa, {\mathcal{M}}, {\mathcal{Y}} \right)
$$

$$
d {\mathcal{X}}/dt = F_{\mathcal{X}} \left(E, {\mathcal{X}} \right)
$$

### Equilibrium

Because the framework is expressed generically, equilibrium solutions can be worked
out for each specific model independently, with input from other components assumed to
be constants. Models can be assembled at equilibrium from these individual model solutions.
Please see the vignettes for each specific model for closed-form solutions.

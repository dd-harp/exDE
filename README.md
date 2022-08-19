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
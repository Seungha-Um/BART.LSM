# Particle MCMC for agent-based SIS model

This repository contains code to implement the methodology described in
the paper “’A Bayesian Approach to Estimate Causal Peer Influence
Accounting for Latent Network Homophil”, by Um, Sweet and Adhikari
(2025)

This package uses the primary functions from
[`SoftBART`](https://github.com/theodds/SoftBART) and the primary
functions from [`HLSM`](https://github.com/cran/HLSM) to incorporate the
SoftBART model and HLSM model as components.

# Installation

The packages can be installed with the `devtools` package:

    library(devtools) 
    devtools::install_github(repo='Seungha-Um/BART.LSM') 

The package with the vignettes can be installed with

    devtools::install_github(repo='Seungha-Um/BART.LSM', build_vignettes = TRUE) 

and then accessed by running `browseVignettes("BART.LSM")`.
Alternatively, vignettes are available at
[Simulation](https://rpubs.com/sheom0808/981711)

Alternatively, vignettes are available at
[Simulation](https://rpubs.com/sheom0808/946709) and [real
data](https://rpubs.com/sheom0808/926959).

One of the vignettes replicates our analysis of the advice-seeking
network dataset, including a noisy subset, while the other demonstrates
our methods using simulated data.

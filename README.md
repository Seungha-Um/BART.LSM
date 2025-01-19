# BART\_LSM model to estimate the average causal peer influence.

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

<!-- and then accessed by running `browseVignettes("BART.LSM")`. Alternatively, vignette is available at [Simulation](https://rpubs.com/sheom0808/981711) -->
<!-- Alternatively,  -->
<!-- vignette is available at [Simulation](https://rpubs.com/sheom0808/946709). -->
<!-- The vignette illustrates our methods using simulated data, while this repository provides the advice-seeking network real dataset, along with a noisy subset. -->

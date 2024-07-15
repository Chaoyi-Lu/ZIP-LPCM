# Simulation Studies for ZIP-LPCM-MFM

This tutorial includes the coding for the two simulation studies illusrated in the paper *A Zero-Inflated Latent Position Cluster Model with Mixture of Finite Mixtures*.

Recall here that we have two different scenarios within each simulation study.
The scenario 1 in the first simulation study focuses on a network randomly generated from a zero-inflated Poisson latent position cluster model (ZIP-LPCM) while the scenario 2 works on a network randomly generated from a Poisson latent positions cluster model (Pois-LPCM).
In our second simulation study, we mainly focus on the networks randomly generated from the zero-inflated Poisson stochastic block model (ZIP-SBM), which is the one we newly proposed in [Lu, C., Durante, D., and Friel, N. [2024+], "Zero-inflated stochastic block modeling of efficiency-security tradeoffs in weighted criminal networks"]().
The synthetic network in scenario 2 is equipped with a hub while the one in scenario 1 is not.

The source code of all the functions required for the experiments in our paper is included in the [`Functions.R`] file.
We load the functions via the code below.

``` r
rm(list=ls())
gc()
source("Functions.R")
```

Within the [`Functions.R`] file, we also added some comments to explain the corresponding lines of the code.
Here we refer to that file for more details.

## Simulation study 1

The scenario 1 network is randomly generated from a ZIP-LPCM with the following code.


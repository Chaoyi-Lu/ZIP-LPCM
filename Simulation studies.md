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

``` r
SS1_Scenario1_Directed_ZIPLPCM <-
  Simulation_Directed_ZIPLPCM(beta=3,P=matrix(c(0.4,0.1,0.05,0.1,0.05,  0.05,0.4,0.1,0.05,0.1,  0.1,0.05,0.4,0.1,0.05,
                                                0.05,0.1,0.05,0.4,0.1,  0.1,0.05,0.1,0.05,0.4),5,5),
                              mu=matrix(c(-1.5,-1.5,-1.5, -2,2,0, 2,-2,0, 2,2,-2, -2,-2,2),nrow=3,ncol=5),
                              tau=c(1/0.25,1/0.5,1/0.75,1/1,1/1.25),d=3,
                              z=c(rep(1,5),rep(2,10),rep(3,15),rep(4,20),rep(5,25)),seed=NULL)
```

The corresponding model parameters and latent variables are in line with those introduced in the paper.
This simulation function brings the following contents for a ZIP-LPCM network.

``` r
SS1_Scenario1_Directed_ZIPLPCM$Y
SS1_Scenario1_Directed_ZIPLPCM$nu
SS1_Scenario1_Directed_ZIPLPCM$z
SS1_Scenario1_Directed_ZIPLPCM$U
SS1_Scenario1_Directed_ZIPLPCM$Z
```

Here, "Y" is the $N \times N$ adjacency matrix of the network, "nu" is the $N \times N$ unusual zero indicator

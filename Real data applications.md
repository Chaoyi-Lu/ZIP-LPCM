# Real data applications for ZIP-LPCM-MFM

This tutorial contains the **zero-inflated Poisson latent position cluster model (ZIP-LPCM)** implementation and post-processing code for all the four **real data applications (RDA)** we illustrate in the paper **"A Zero-Inflated Latent Position Cluster Model with Mixture of Finite Mixtures" (ZIP-LPCM-MFM)**.
The functions required during the experiments are all included in the [`Functions.R`] file of this repository and can be loaded via:

``` r
rm(list=ls())
gc()
source("Functions.R")
```

The 1st real network is the **Sampson Monks** network which is publicly available in the `latentnet` package following:

``` r
# Sampson Monks network, Directed
library(latentnet)
data(sampson)
```

The adjacency matrix of this network can be extracted by:

``` r
SampsonMonks_Directed_adj <- as.matrix(samplike,attrname='nominations') # Nominations adjacency matrix
colnames(SampsonMonks_Directed_adj ) <- c()
rownames(SampsonMonks_Directed_adj ) <- c()
```

The reference clustering provided by Sampson contains three groups: "Turks", "Outcasts" and "Loyal". 
The clustering is available and can be transformed to numbers following:

``` r
sampson_monks_group <- as.integer(c(as.factor(samplike%v%"group")))
```

Recall here that we propose to apply the [Rastelli, R. and Friel, N. (2018)](https://pubmed.ncbi.nlm.nih.gov/30220822/) label-switching method for all the clustering, i.e., assigning the node 1 to group 1 by default, and then iteratively assigning the next node either to a new empty group or an existing group.
Thus the **Sampson Monks** reference grouping is label-switched following:

``` r
sampson_monks_group <- LabelSwitching_SG2003(z=sampson_monks_group)$z
```

Each number $1,2,3$ in the the reference clustering is one-to-one correspondance to the "Turks", "Outcasts" and "Loyal" group of the **Sampson Monks** shown as:

``` r
sampson_monks_group
samplike%v%"group"
```

Recall also here that we treat the combination of the reference grouping and an indicator of whether or not each monk attended the minor seminary of "Cloisterville" before coming to the monastery as the exogenous node attributes leveraged to apply supervised implementations in practice.
The "Cloisterville" information is available via:

``` r
samplike%v%"cloisterville" # An indicator of attendance the minor seminary of "Cloisterville" before coming to the monastery.
```

And we combine such information with the reference grouping for each monk following:

``` r
sampson_monks_group_cloisterville <- paste(samplike%v%"group",samplike%v%"cloisterville",sep="_")
# sampson_monks_group_cloisterville
sampson_monks_group_cloisterville_NodeA <- as.integer(c(as.factor(sampson_monks_group_cloisterville)))
sampson_monks_group_cloisterville_NodeA_temp <- sampson_monks_group_cloisterville_NodeA
sampson_monks_group_cloisterville_NodeA[sampson_monks_group_cloisterville_NodeA_temp==6] <- 1
sampson_monks_group_cloisterville_NodeA[sampson_monks_group_cloisterville_NodeA_temp==4] <- 3
sampson_monks_group_cloisterville_NodeA[sampson_monks_group_cloisterville_NodeA_temp==2] <- 5
sampson_monks_group_cloisterville_NodeA[sampson_monks_group_cloisterville_NodeA_temp==5] <- 2
sampson_monks_group_cloisterville_NodeA[sampson_monks_group_cloisterville_NodeA_temp==3] <- 4
sampson_monks_group_cloisterville_NodeA[sampson_monks_group_cloisterville_NodeA_temp==1] <- 6
sampson_monks_group_cloisterville_NodeA
```

which brings our exogenous node attributes we used in practice.


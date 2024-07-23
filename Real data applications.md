# Real data applications for ZIP-LPCM-MFM

This tutorial contains the **zero-inflated Poisson latent position cluster model (ZIP-LPCM)** implementation and post-processing code for all the four **real data applications (RDA)** we illustrate in the paper **"A Zero-Inflated Latent Position Cluster Model with Mixture of Finite Mixtures" (ZIP-LPCM-MFM)**.
The functions required during the experiments are all included in the [`Functions.R`] file of this repository and can be loaded via:

``` r
rm(list=ls())
gc()
source("Functions.R")
```

## 1. Loading real networks

### 1.1 Sampson monks

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

The [Rastelli, R. and Friel, N. (2018)](https://pubmed.ncbi.nlm.nih.gov/30220822/) label-switching method is applied for the clustering, i.e., assigning the node 1 to group 1 by default, and then iteratively assigning the next node either to a new empty group or an existing group.
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

### 1.2 Windsurfers and Train bombing

Both of the **Windsurfers** and the **Train bombing** real networks are publicly avaiable at the website [http://konect.cc/networks/](http://konect.cc/networks/).
We also upload the real data at [`Datasets/windsurfers.txt`] and [`Datasets/train_bombing.txt`] of this repository, and the networks can be directly extracted by:

``` r
Windsurfers <- read.table("Datasets/windsurfers.txt", skip=2, sep="")
Train_bombing <- read.table("Datasets/train_bombing.txt", skip=2, sep="")
```

The corresponding adjacency matrices can be obtained following:

``` r
# Windsurfers Network Dataset, UnDirected
library(igraph) 
Windsurfers_graph <- graph_from_edgelist(as.matrix(Windsurfers[,1:2]))
Windsurfers_graph <- set_edge_attr(graph=Windsurfers_graph, name="Weight", value=c(Windsurfers[,3]))
Windsurfers_adj <- as.matrix(get.adjacency(Windsurfers_graph,attr = "Weight"))
library(gdata)
lowerTriangle(Windsurfers_adj,byrow=TRUE) <- upperTriangle(Windsurfers_adj)

# Train bombing Network Dataset, UnDirected
library(igraph) 
Train_bombing_graph <- graph_from_edgelist(as.matrix(Train_bombing[,1:2]))
Train_bombing_graph <- set_edge_attr(graph=Train_bombing_graph, name="Weight", value=c(Train_bombing[,3]))
Train_bombing_adj <- as.matrix(get.adjacency(Train_bombing_graph,attr = "Weight"))
library(gdata)
lowerTriangle(Train_bombing_adj,byrow=TRUE) <- upperTriangle(Train_bombing_adj)
```

### 1.3 The 'Ndrangheta Mafia Network

The original data is available at [https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndrangheta-mafia-2](https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndrangheta-mafia-2) and we follow the same pre-processing steps as [Legramanti et al. (2022)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-4/Extended-stochastic-block-models-with-application-to-criminal-networks/10.1214/21-AOAS1595.short) to process the raw data and to focus on the processed data for implementations.
We refer to the [GitHub page](https://github.com/danieledurante/ESBM/blob/master/Application/application.md) of the paper Legramanti et al. (2022) for more details, and the following pre-processing code are extracted from their GitHub page with some minor modifications on variable names.

``` r
A <- read.csv(file="Datasets/NDRANGHETAMAFIA_2M.csv",header=TRUE, stringsAsFactors = TRUE)
A[23,20] <- 1
A[48,23] <- 1
A[78,35] <- 1
A[115,22] <- 1
# Suspects who never attended a summit
sel_empty <- which(apply(as.matrix(A[,-1]),1,sum)==0)
# Suspects who have not been recognized during the investigation process
A[c(38,105,106,125,135),1]
# Indicators of the two groups of suspects to be excluded
sel_empty <- c(c(sel_empty),c(38,105,106,125,135))
# remove these suspects from the dataset
A <- A[-sel_empty,]
# Create a vector with the actors' names
actors <- A[,1]
actors <- droplevels(actors)
A <- as.matrix(A[,-1])
A <- A%*%t(A)
diag(A) <- 0
dim(A)
rownames(A) <- colnames(A) <- c(1:dim(A)[1])
# ----------------------------------------------------
# Locale membership ("OUT": Suspects not belonging to La Lombardia. "MISS": Information not available)
Locale <- c("C","OUT","A","MISS","O","A","MISS","D","D","D","D","D","C","P","L","L","Q","MISS","B","OUT",
            "B","B","I","MISS","OUT","D","A","O","N","N","H","OUT","D","E","G","G","L","A","OUT","Q",
            "C","OUT","Q","L","C","MISS","C","C","F","C","OUT","D","A","B","B","E","M","MISS","C","C",
            "C","B","H","C","C","E","E","E","E","C","MISS","L","A","A","E","E","C","E","E","E",
            "C","MISS","OUT","C","C","E","G","A","A","B","I","I","A","B","B","OUT","I","A","G","N",
            "E","D","F","OUT","OUT","C","D","C","MISS","MISS","C","MISS","E","E","C","MISS","OUT","B","L","A",
            "D","D","O","MISS","B","D","O","D","D","A","A","I","C","MISS","MISS","MISS","A","A","F","E",
            "C","Q","H","B","B","B")

# ----------------------------------------------------
# Leadership role ("miss": Information not available)
Role <- c("aff","aff","aff","miss","aff","boss","miss","boss","boss","aff","aff","aff","aff","aff","aff","aff","aff","miss","aff","boss",
          "boss","boss","boss","miss","boss","boss","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","aff","aff","boss","aff",
          "aff","aff","boss","aff","aff","miss","aff","aff","aff","aff","boss","aff","aff","aff","aff","aff","boss","miss","aff","aff",
          "aff","aff","boss","boss","boss","aff","boss","aff","aff","boss","miss","aff","aff","boss","boss","aff","aff","aff","aff","aff",
          "aff","miss","boss","aff","aff","aff","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","boss","boss","boss","boss",
          "aff","aff","aff","aff","aff","aff","aff","boss","miss","miss","aff","miss","aff","aff","aff","miss","aff","aff","boss","aff",
          "aff","aff","aff","miss","aff","aff","boss","aff","boss","aff","aff","aff","aff","miss","miss","miss","aff","aff","boss","aff",
          "aff","aff","boss","aff","aff","boss")

# Indicators of those suspects who are not known to be part of the organization La Lombardia
sel_miss <- which(Locale=="MISS" | Locale=="OUT")
# Remove such suspects from the dataset
Locale_temp <- Locale[-sel_miss]
Role_temp <- Role[-sel_miss]
actors_temp <- actors[-sel_miss]
actors_temp <- droplevels(actors_temp)
A <- A[-sel_miss,-sel_miss]
rownames(A) <- colnames(A) <- c(1:dim(A)[1])
dim(A)
# Interaction between Locale and Role
RoleLocale_temp <- paste(Role_temp,Locale_temp,sep="_")

# Create the dataset used for modeling and inference
sel <- which(Locale_temp=="A" | Locale_temp=="B" | Locale_temp=="C" | Locale_temp=="D" | Locale_temp=="E")
CriminalNetwork <- A[sel,sel]
RoleLocale <- RoleLocale_temp[sel]
Role <- Role_temp[sel]
Locale <- Locale_temp[sel]
actors <- actors_temp[sel]
actors <- droplevels(actors)
rownames(CriminalNetwork) <- colnames(CriminalNetwork) <- c(1:dim(CriminalNetwork)[1])

# Define the node attributes
RDA_criminalNet_A <- c(as.factor(RoleLocale))
RDA_criminalNet_A<- as.integer(RDA_criminalNet_A)
RDA_criminalNet_A[which(RDA_criminalNet_A>5)] <- 7
# Actors with known peripheral roles in locale D
RDA_criminalNet_A[c(6,8,69,73)]<- 6
```

The `Rolelocale` above corresponds to the reference clustering throughout the experiments.
We store the above processed data by:

``` r
write.csv(CriminalNetwork,"Datasets/CriminalNet_Y.csv", row.names = FALSE)
write.csv(RoleLocale,"Datasets/CriminalNet_RoleLocale.csv", row.names = FALSE)
write.csv(actors,"Datasets/CriminalNet_actors.csv", row.names = FALSE)
write.csv(Locale,"Datasets/CriminalNet_Locale.csv", row.names = FALSE)
write.csv(Role,"Datasets/CriminalNet_Role.csv", row.names = FALSE) # "affiliate" or "Boss"
write.csv(RDA_criminalNet_A,"Datasets/CriminalNet_A.csv", row.names = FALSE)
```

The data is uploaded at [`Datasets/`] file of this repository and can be directed loaded via:

``` r
# 'Ndrangheta Mafia Network, Undirected
RDA_criminalNet <- list(Y = as.matrix(read.csv("Datasets/CriminalNet_Y.csv",header = TRUE)),
                        RoleLocale = c(as.matrix(read.csv("Datasets/CriminalNet_RoleLocale.csv",header = TRUE))),
                        Actors = c(as.matrix(read.csv("Datasets/CriminalNet_actors.csv",header = TRUE))),
                        Locale = c(as.matrix(read.csv("Datasets/CriminalNet_Locale.csv",header = TRUE))),
                        role = c(as.matrix(read.csv("Datasets/CriminalNet_Role.csv",header = TRUE))),
                        A = c(as.matrix(read.csv("Datasets/CriminalNet_A.csv",header = TRUE))))
colnames(RDA_criminalNet$Y) <- c()
```

## 2. Implementations and post processing

For each real network, the inference algorithm was implemented for 60,000 iterations with 30,000-iteration burn-in in order for sufficient mixing.
All the **RDA** output shown in the **ZIP-LPCM-MFM** paper are reproducible by setting the random number generator (RNG) seed by `set.seed(1)`.
Multiple implementations can be easily applied by setting different seeds for different rounds of implementations or by simply removing the seed via `set.seed(NULL)` if it's not required to reproduce the output.

### 2.1 Sampson Monks








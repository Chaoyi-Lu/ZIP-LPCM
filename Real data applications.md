# Real Data Applications for ZIP-LPCM-MFM

This tutorial contains the **zero-inflated Poisson latent position cluster model (ZIP-LPCM)** implementation and post-processing code for all the four **real data applications (RDA)** we illustrate in the paper **"A Zero-Inflated Latent Position Cluster Model with Mixture of Finite Mixtures" (ZIP-LPCM-MFM)**.
The functions required during the experiments are all included in the [`Functions.R`] file of this repository and can be loaded via:

``` r
rm(list=ls())
gc()
source("Functions.R")
```

## 1. Loading Real Networks

### 1.1 Sampson Monks Network

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

### 1.2 Windsurfers and Train Bombing Networks

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

### 1.3 'Ndrangheta Mafia Network

The original data is available at [https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndrangheta-mafia-2](https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndrangheta-mafia-2) and we follow the same pre-processing steps as [Legramanti et al. (2022)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-4/Extended-stochastic-block-models-with-application-to-criminal-networks/10.1214/21-AOAS1595.short) to process the raw data and to focus on the processed data for implementations.
We refer to the [GitHub page](https://github.com/danieledurante/ESBM/blob/master/Application/application.md) of the paper Legramanti et al. (2022) for more details, and the following pre-processing code are extracted from their GitHub page with some minor modifications on variable names.

``` r
A <- read.csv(file="Datasets/NDRANGHETAMAFIA_2M.csv",header=TRUE, stringsAsFactors = TRUE)
# Correct the recorded data based on the judicial acts
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

## 2. Implementations and Post-processing

For each real network, the inference algorithm was implemented for 60,000 iterations with 30,000-iteration burn-in in order for sufficient mixing.
All the **RDA** output shown in the **ZIP-LPCM-MFM** paper are reproducible by setting the seed of the random number generator (RNG) as `set.seed(1)`.
Multiple implementations can be easily applied by setting different seeds for different rounds of implementations or by simply removing the seed via `set.seed(NULL)` if it's not required to reproduce the output.

### 2.1 Sampson Monks Network

The **Sampson Monks** network is a directed network and the presumed exogenous node attributes are avaiable in **section 1.1** above.
So **supervised ZIP-LPCM** is implemented in practice for such a real network following:

``` r
# Sampson monks directed real network supervised ZIP-LPCM T = 60000 round 1
set.seed(1)
start.time <- Sys.time()
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1 <- 
  MwG_Directed_ZIPLPCM(Y = SampsonMonks_Directed_adj,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=9,
                       sigma2prop_beta=0.4^2,sigma2prop_U=0.2,d=3,z=1:nrow(SampsonMonks_Directed_adj),
                       p_eject=0.5,A=sampson_monks_group_cloisterville_NodeA,omega_c=1)
end.time <- Sys.time()
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_time <- end.time - start.time
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_time
# Time difference of 34.79106 mins
```

where a reference running time is provided above for a laptop with with eight 1.80GHz processors.
The prior parameters are similar to those implemented in the simulation studies, and the proposal variances are tuned so that the acceptance rates of the Metropolis-Hastings steps of $\boldsymbol{U}$ and $\beta$ are:

``` r
# Define the burn in
iteration_after_burn_in <- 30002:60001
burn_in <- 30001

## Check U acceptance rate
apply(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)
## Check the mean of U acceptance rate
mean(apply(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)) # 0.18595
## Check beta acceptance rate
mean(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$acceptance_count_beta[iteration_after_burn_in-1]) # 0.2164333
```

The post-processing steps begin with applying label-switching on the posterior clustering and posterior unusual zero probability:

``` r
# Apply label switching on the post clustering z and post unusual zero probability P
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz <- matrix(NA,nrow=nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z),ncol=ncol(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z))
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSP <- list()
for (t in 1:nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z)){
  LS_temp <- LabelSwitching_SG2003(z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z[t,],
                                   matrix_KbyK = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$P[[t]])
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[t,] <- LS_temp$z
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSP[[t]] <- LS_temp$matrix_KbyK
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
```

and with applying Procrustes transform on posterior latent positions with respect to the last posterior sample of the latent positions:

``` r
# Apply Procrustes Transform on posterior latent positions U
library("IMIFA")
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_PTU <- list()
for (t in 1:nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z)){
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_PTU[[t]] <-
    Procrustes(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$U[[t]],
               RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$U[[60001]], translate = TRUE ,dilate = FALSE)$X.new
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
```

Then we can first check the sufficient mixing of posterior samples via the traceplot of the complete likelihoods and the number of clusters:

``` r
## Check complete likelihood for each iteration
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_Like <- c()
for (t in 1:nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z)){
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_Like <-
    c(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_Like,
      Directed_ZIPLPCM_MFM_CompleteLikelihood(X=RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$X[[t]],
                                              U=RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_PTU[[t]],
                                              beta=RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$beta[t],
                                              nu=RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$nu[[t]],
                                              P=RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSP[[t]],
                                              z=RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[t,],
                                              alpha1=1,alpha2=0.103,omega=0.01,  alpha=3,
                                              A=sampson_monks_group_cloisterville_NodeA,omega_c=1))
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
plot(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_Like,type = "l",xlab = "",ylab = "", main = "Likelihood",cex.axis = 0.8)

# Check the trace plot of K
plot(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$K,type = "l",xlab = "",ylab = "", main = "K Trace Plot",cex.axis = 0.8)
```

To obtain a point estimate of the posterior clustering, we follow equation (29) of the **ZIP-LPCM-MFM** paper and leverage the greedy algorithm proposed by [Rastelli, R. and Friel, N. (2018)](https://pubmed.ncbi.nlm.nih.gov/30220822/).
However, we start from obtaining the marginal posterior mode of the posterior clustering and treat it as the initial state of the Rastelli, R. and Friel, N. (2018) greedy method:

``` r
# Summarize posterior clustering z by the greedy algorithm proposed by Rastelli and Friel (2018)
# We start from obtaining the marginal posterior mode of the posterior z chain
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States <- c() # initialize a list which will store different clustering states; the clustering states are labeled from 1,2,3... and are put at the 1st,2nd,3rd... row of the matrix, respectively
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration <- list() # store all the t's (iteration number) which provides the same cluster as the clustering state 1,2,3...
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop <- iteration_after_burn_in # the "IterationLoop" which stores all the iteration t's which we focus on, that is, all the iteration t's after burn-in
StatesLabelIndicator = 0 # initialize the label for the clustering states
while (length(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop)!=0){ # if the "IterationLoop" is not empty
  StatesLabelIndicator <- StatesLabelIndicator + 1 # assign the next label to the next clustering state
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_FirstState <- RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop[1],] # extract the first clustering state for the "IterationLoop"
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States <- rbind(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States,
                                                                             RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_FirstState) # store the first state within the "IterationLoop" with label "StatesLabelIndicator" in the list which will contain all different unique states
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp <- c() # create a vector to temporarily store all the iteration t's whose clustering is the same as the first clustering state within the "IterationLoop"
  for (t in RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop){ # loop over all the current existing iterations in "IterationLoop"
    if (sum(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[t,]==RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_FirstState)==length(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_FirstState)){ # if the t's clustering is the same as the "FirstState"
      RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp <- c(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp,t) # store the iteration t in the temporary vector
    }
  }
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration[[StatesLabelIndicator]] <- RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp # store all the t's as the list element
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop <-
    (iteration_after_burn_in)[-(unlist(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration)-burn_in)] # remove all the iterations we have stored and then move to the next clustering state
}
rownames(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States) <- NULL
nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States) # check the number of different clustering states
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency <- c() # check the number of times each clustering state occurs
for (t in 1:nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States)){
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency <-
    c(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency,
      length(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration[[t]]))
}
library(GreedyEPL) # obtain the summarized z by the greedy algorithm proposed by Rastelli and Friel (2018)
output <- MinimiseEPL(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,],
                      list(Kup = 20, loss_type = "VI",decision_init = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States[
                        which.max(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency),]))
output$EPL # check minEVI posterior loss: 0.03696773
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z <- output$decision # Save output as summarized z
table(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z,sampson_monks_group,dnn = c("","")) # compare to the reference clustering
```

It's shown that the summarized clustering is exactly the same as the reference clustering.
Based on the results we obtained above, we can also follow equation (30) of the **ZIP-LPCM-MFM** paper to obtain a point estimate $\hat{\boldsymbol{U}}$ for the latent positions:

``` r
# Obtain summarized U
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_hat_zIteration <- # extract all the iterations whose clustering is identical to hat_z
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration[[
    which(apply(t(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_States)==
                  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z,2,sum)==length(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z))
  ]]
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_MaxLikeHatz_iteration <- RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_hat_zIteration[
  which.max(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_Like[RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz_hat_zIteration])]
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U <-
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_PTU[[RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_MaxLikeHatz_iteration]]
```

Once we obtained both $\hat{\boldsymbol{z}}$ and $\hat{\boldsymbol{U}}$, the 3-dimensional interactive plot of them can be produced following:

``` r
library("igraph")
library("RColorBrewer")
library("ggplot2")
My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])
g_obs <- graph_from_adjacency_matrix(SampsonMonks_Directed_adj,mode = "directed",weighted = TRUE)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(SampsonMonks_Directed_adj))[E(g_obs)$weight]
betw <- betweenness(g_obs) # evaluate the betweeness of the network
VertexSize <- sqrt(betw/1.5+mean(betw))*1 # set the vertex size
library("plotly")
fig <- plot_ly() %>% # plot the summarized hat_z and hat_U
  add_markers(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(SampsonMonks_Directed_adj),"<br>z*:",sampson_monks_group),
              size=VertexSize,sizes=c(300,600),
              color=as.factor(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z),colors=My_colors[6:8]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.75*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "hat_U and hat_z",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

which is uploaded at [`/Interactive 3-d latent positions plots/RDA_SampsonMonks_InteractivePlot.html`] of this repository.
Since the GitHub page cannot directly view this file, we suggest the readers to download it first and then the 3-d interactive plot would be free to play with.
Note that running the code above would bring some warning messages and these are raised by some package bugs which are not yet resolved.
However, these bugs would not affect our output and the readers can ignore those warning messages.
The above code for the interactive 3-d plot also attached some texts or notes for each node and each non-zero interaction to help readers have a better understanding of the network.
If the readers put the mouse pointer on each node of the interactive plot, there will be a comment bracket showing (i) the coordinate of the node, (ii) the node number (e.g. node 1, node 2, ...), (iii) the reference clustering of the node.
If the mouse pointer is put on each non-zero interaction, the bracket will show (i) either the start coordinate or the end coordinate of the interaction vector, (ii) the interaction weight, (iii) an indicator of whether the interaction is from node $i$ to node $j$ where $i\< j$, i.e., whether the corresponding $y_{ij}$ is the upper-diagonal entry of the $\boldsymbol{Y}$ for the directed networks.

The **Figure 7** of the **ZIP-LPCM-MFM** paper can thus be reproduced by the `subplot()` and the `orca()` functions following the code below.
The `subplot()` function helps us add two interactive plots with different default viewing angle into one figure while the `orca()` function can produce a high quality screenshot of the interactive plots.
Directly taking a screenshot of the interactive plots would bring a low-quality figure which is not satisfactory.
However, the readers need to first download the "orca" app following [https://plotly.com/r/static-image-export/](https://plotly.com/r/static-image-export/) and [https://github.com/plotly/orca#installation](https://github.com/plotly/orca#installation) in order to use the `orca()` function in the `plotly` package.

``` r
# Use subplot() to plot two interactive plot with different viewing angle
library("plotly")
# Plot the 1st viewing angle of the latent positions
fig1 <- plot_ly(scene ="scene1") %>% # plot the summarized clustering and hat_U
  add_markers(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[,1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[,2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(SampsonMonks_Directed_adj),"<br>c:",sampson_monks_group_cloisterville),
              size=VertexSize,sizes=c(300,600),showlegend = FALSE,
              color=as.factor(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_z),colors=My_colors[6:8]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig1 <- fig1 %>%
    add_trace(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.75*E(g_obs)$weight[i]))
}


# Plot the 2nd viewing angle of the latent positions
fig2 <- plot_ly(scene ="scene2") %>% # plot the summarized clustering and hat_U
  add_markers(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[,1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[,2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(SampsonMonks_Directed_adj),"<br>c:",sampson_monks_group_cloisterville),
              size=VertexSize,sizes=c(300,600),showlegend = FALSE,
              color=as.factor(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_z),colors=My_colors[6:8]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig2 <- fig2 %>%
    add_trace(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.75*E(g_obs)$weight[i]))
}

Fig <- subplot(fig1, fig2)
Fig <- Fig %>% layout(title = "", margin = list(l = 0,r = 0,b = 0,t = 0,pad = 0),
                        scene = list(domain=list(x=c(0,1/2),y=c(0,1)),
                                     xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                     camera = list(eye = list(x = 0.1, y = -1.75, z = 0.75)),
                                     aspectmode='auto'),
                        scene2 = list(domain=list(x=c(1/2,0.999),y=c(0,1)),
                                      xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                      camera = list(eye = list(x = 1.75, y = 0.1, z = 0.75)),
                                      aspectmode='auto'))
Fig # directly make the plot in Rstudio
# Apply orca() to obtain a high quality screenshot
orca(Fig, "RDA_SampsonMonks_hat_U_Y.pdf",scale=1,width=1800,height=850)
```

The summarized $\hat{\beta}$ obtained by posterior mean can also be easily obtained:

``` r
## Summarize beta
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_beta <-
  mean(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$beta[iteration_after_burn_in])
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_beta #  1.008339
```

As no corresponding reference to compare with, we propose not to show the $\hat{\beta}$ in the paper.
But we provide above the $\hat{\beta}$ value we obtained.

### 2.2 Windsurfers Network

### 2.3 Train Bombing Network

### 2.4 'Ndrangheta Mafia Network





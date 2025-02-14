# Real Data Applications for ZIP-LPCM

This tutorial contains the **zero-inflated Poisson latent position cluster model (ZIP-LPCM)** implementation and post-processing code for all the four **real data applications (RDA)** we illustrate in the paper **"A Zero-Inflated Poisson Latent Position Cluster Model" (ZIP-LPCM)**.
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

### 1.3 Summit co-attendance criminality Network

The original data was processed by [Legramanti et al. (2022)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-4/Extended-stochastic-block-models-with-application-to-criminal-networks/10.1214/21-AOAS1595.short), and we focus on the processed data available at [GitHub page](https://github.com/danieledurante/ESBM/blob/master/Application/application.md) for implementations.
Following the processing code therein with just changing the variable names, the processed data is also available at [`Datasets/`] file of this repository and can be directly loaded via:

``` r
# Summit co-attendance criminality Network, Undirected
RDA_criminalNet <- list(Y = as.matrix(read.csv("Datasets/CriminalNet_Y.csv",header = TRUE)),
                        RoleAffiliation = c(as.matrix(read.csv("Datasets/CriminalNet_RoleAffiliation.csv",header = TRUE))),
                        Affiliation = c(as.matrix(read.csv("Datasets/CriminalNet_Affiliation.csv",header = TRUE))),
                        role = c(as.matrix(read.csv("Datasets/CriminalNet_Role.csv",header = TRUE))),
                        A = c(as.matrix(read.csv("Datasets/CriminalNet_A.csv",header = TRUE))))
colnames(RDA_criminalNet$Y) <- c()
```

## 2. Implementations and Post-processing

For each real network, the inference algorithm was implemented for 60,000 iterations with 30,000-iteration burn-in in order for sufficient mixing.
All the **RDA** output shown in the **ZIP-LPCM** paper are reproducible by setting the seed of the random number generator (RNG) as `set.seed(1)`.
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

To obtain a point estimate of the posterior clustering, we follow equation (29) of the **ZIP-LPCM** paper and leverage the greedy algorithm proposed by [Rastelli, R. and Friel, N. (2018)](https://pubmed.ncbi.nlm.nih.gov/30220822/).
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
Based on the results we obtained above, we can also follow equation (30) of the **ZIP-LPCM** paper to obtain a point estimate $\hat{\boldsymbol{U}}$ for the latent positions:

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
              stroke = I("black"), span = I(1),
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
Since in this real data application the inferred $\hat{\boldsymbol{z}}$ is exactly the same as the reference clustering, so we replace the reference clustering shown in such comment brackets with the node attributes `sampson_monks_group_cloisterville` in the rest 3-d plots shown next.
If the mouse pointer is put on each non-zero interaction, the bracket will show (i) either the start coordinate or the end coordinate of the interaction vector, (ii) the interaction weight, (iii) an indicator of whether the interaction is from node $i$ to node $j$ where $i\< j$, i.e., whether the corresponding $y_{ij}$ is the upper-diagonal entry of the $\boldsymbol{Y}$ for the directed networks.

The **Figure 8** of the **ZIP-LPCM** paper can thus be reproduced by the `subplot()` and the `orca()` functions following the code below.
The `subplot()` function helps us add two interactive plots with different default viewing angle into one figure while the `orca()` function can produce a high quality screenshot of the interactive plots.
Directly taking a screenshot of the interactive plots would bring a low-quality figure which is not satisfactory.
However, the readers need to first download the "orca" app following [https://plotly.com/r/static-image-export/](https://plotly.com/r/static-image-export/) and [https://github.com/plotly/orca#installation](https://github.com/plotly/orca#installation) in order to use the `orca()` function in the `plotly` package.

``` r
# Use subplot() to plot two interactive plot with different viewing angle
library("plotly")
# Plot the 1st viewing angle of the latent positions
fig1 <- plot_ly(scene ="scene1") %>% # plot the summarized clustering and hat_U
  add_markers(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(SampsonMonks_Directed_adj),"<br>c:",sampson_monks_group_cloisterville),
              size=VertexSize,sizes=c(300,600),showlegend = FALSE,
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z),colors=My_colors[6:8]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig1 <- fig1 %>%
    add_trace(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.75*E(g_obs)$weight[i]))
}


# Plot the 2nd viewing angle of the latent positions
fig2 <- plot_ly(scene ="scene2") %>% # plot the summarized clustering and hat_U
  add_markers(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(SampsonMonks_Directed_adj),"<br>c:",sampson_monks_group_cloisterville),
              size=VertexSize,sizes=c(300,600),showlegend = FALSE,
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z),colors=My_colors[6:8]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig2 <- fig2 %>%
    add_trace(x = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
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

Recall here that the posterior mean of the unusual zero indicator $\boldsymbol{\nu}$ approximates the conditional probability of unusual zeros provided that the corresponding observed interactions are zeros, i.e., the equation (22) of the **ZIP-LPCM** paper.
We can obtain such a statistic by:

``` r
## Obtain posterior mean of nu, i.e. approximate P_m0
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_nu <-
  Reduce("+",RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$nu[iteration_after_burn_in])/length(iteration_after_burn_in)
```

Here we use `P_m0` in the code to denote such a conditional probability.

Similar to the simulation studies shown in the [`Simulation studies.md`] file of this repository, we can also obtain for each posterior state the individual-level unusual zero probability $\boldsymbol{p}$ which is a $N \times N$ matrix with each entry $i,j$ being $p_{z_iz_j}$ for $i,j = 1,2,\dots,N$ following the code below.
And we can thus obtain the posterior mean of the posterior samples of $\boldsymbol{p}$ to obtain the corresponding summary statistic accounting for the uncertainty of the posterior clustering.

``` r
# Obtain the individual-level unusual zero probability p for each iteration
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_p <- list()
library(Rfast) # for Rfast::Dist()
for (t in 1:nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[t,],length(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[t,]),
                  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$K[t]))==(1:RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$K[t]))*1
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_p[[t]] <- Z%*%RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSP[[t]]%*%t(Z)
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p <- Reduce("+",RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_p[iteration_after_burn_in])/length(iteration_after_burn_in)
diag(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p) <- 0
```

Similarly, the individual-level Poisson rate matrix with each entry being $\lambda_{ij}=\text{exp}(\beta-||\boldsymbol{u_i}-\boldsymbol{u_j}||)$ can also be obtained for each posterior state, and the corresponding posterior mean summary statistic can thus be obtained:

``` r
# Obtain the individual-level Poisson rate lambda for each iteration
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_lambda <- list()
library(Rfast) # for Rfast::Dist()
for (t in 1:nrow(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[t,],length(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_LSz[t,]),RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$K[t]))==(1:RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$K[t]))*1
  Dist_U <- Rfast::Dist(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_PTU[[t]])
  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_lambda[[t]] <- exp(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1$beta[t]-Dist_U)
  if ((t%%1000) == 0){cat("t=",t,"\n") }
}
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_lambda <- Reduce("+",RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_lambda[iteration_after_burn_in])/length(iteration_after_burn_in)
diag(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_lambda) <- 0
```

Based on the above individual-level statistics, the heatmap plots shown in **Figure 7** of the **ZIP-LPCM** paper can be recovered by:

``` r
library("RColorBrewer")
library("pheatmap")
library("ggplot2")

My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

SampsonMonks_Directed_adj_dataframe <- as.data.frame(SampsonMonks_Directed_adj)
rownames(SampsonMonks_Directed_adj_dataframe) <- colnames(SampsonMonks_Directed_adj_dataframe) <- 1:nrow(SampsonMonks_Directed_adj)

RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe <- as.data.frame(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_nu)
rownames(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe) <- colnames(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe) <- 1:nrow(SampsonMonks_Directed_adj)

# Plot of the adj, adj|hat_z and hat_nu*(1-dpois(0,hat_lambda))|hat_z
annotation_row_z_ref <- as.data.frame(as.factor(matrix(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z,length(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z),1)))
colnames(annotation_row_z_ref) <- "z_ref"
annotation_colors_z_ref <- My_colors[6:8]
names(annotation_colors_z_ref) <- sort(unique(annotation_row_z_ref$z_ref))
annotation_colors_blank <- rep(brewer.pal(9,"Greys")[1],3)
names(annotation_colors_blank) <- sort(unique(annotation_row_z_ref$z_ref))

RDA_SampsonMonks_Directed_Y_heatmap <-
  pheatmap(SampsonMonks_Directed_adj_dataframe,
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(max(SampsonMonks_Directed_adj))),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_blank),
           annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE)

RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_Y_hat_z_heatmap <-
  pheatmap(SampsonMonks_Directed_adj_dataframe[order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z),order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(max(SampsonMonks_Directed_adj))),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)))

RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_P_m0_Non0_heatmap <-
  pheatmap((RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe*(1-dpois(0,RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_lambda)))[order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z),order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:4])(100)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)))

library("grid")
library("gridExtra")
g <- grid.arrange(RDA_SampsonMonks_Directed_Y_heatmap[[4]],
                  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_Y_hat_z_heatmap[[4]],
                  RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_P_m0_Non0_heatmap[[4]],
                  nrow=1,ncol=3,vp=viewport(width=1, height=1))
Fig <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(colour = "white", linewidth = 0))
Fig
```

And the heatmap plots of the individual-level unusual zero probability $\hat{\boldsymbol{p}}$ and the approximate `P_m0`, i.e., the $\hat{\boldsymbol{\nu}}$ which are not illustrated in the paper can also be obtained here:

``` r
# Heatmap plot of the hat_nu which approximates the P_m0
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_P_m0_heatmap <-
  pheatmap(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe[order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z),order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:5])(100)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)))
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_P_m0_heatmap_draw <- cowplot::ggdraw(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_P_m0_heatmap[[4]])+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[1]))
print(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_P_m0_heatmap_draw)

# Heatmap plot of the individual-level hat_p
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe <- as.data.frame(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p)
rownames(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe) <- colnames(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe) <- 1:nrow(SampsonMonks_Directed_adj)

RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap <-
  pheatmap(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe[order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_z),order(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:3])(100)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_R1_hat_z))!=0)))
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap_draw <- cowplot::ggdraw(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap[[4]])+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[1]))
print(RDA_SampsonMonks_Directed_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap_draw)
```

### 2.2 Windsurfers Network

We follow similar steps as Section 2.1 above to analyze the **Windsurfers** real network.
This real network is a undirected network and no reference clustering or exogenous node attributes is available to us.
Thus an unsupervised implementation of the ZIP-LPCM analysis is applied in practice.
Further, since the size of this network is significantly larger than the **Sampson Monks** network, we propose $\text{Beta}(1,9)$ unusual zero probability prior to encourage the clustering.
Once again, all the output are reproducible by setting `set.seed(1)`:

``` r
# Windsurfers undirected real network unsupervised ZIP-LPCM T = 60000 round 1
set.seed(1)
start.time <- Sys.time()
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1 <- 
  MwG_UnDirected_ZIPLPCM(Y = Windsurfers_adj,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=19,
                         sigma2prop_beta=0.15^2,sigma2prop_U=0.2,d=3,z=1:nrow(Windsurfers_adj),
                         p_eject=0.5)
end.time <- Sys.time()
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_time <- end.time - start.time
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_time
# Time difference of 1.2212 hours
```

It's shown that larger network size brings higher computational burden, and the implementation in our experiment took 1.2212 hours to finish the running.

The following post-processing steps are also similar to the **Sampson Monks** network case shown in Section 2.1:

``` r
# Define the burn in
iteration_after_burn_in <- 30002:60001
burn_in <- 30001

## Check U acceptance rate
apply(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)
mean(apply(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)) # 0.236455
## Check beta acceptance rate
mean(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$acceptance_count_beta[iteration_after_burn_in-1]) # 0.2179667

# Apply label switching on the post clustering z
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz <- matrix(NA,nrow=nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z),ncol=ncol(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z))
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP <- list()
for (t in 1:nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  LS_temp <- LabelSwitching_SG2003(z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z[t,],
                                   matrix_KbyK = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$P[[t]])
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,] <- LS_temp$z
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP[[t]] <- LS_temp$matrix_KbyK
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
# Apply Procrustes Transform on posterior latent positions U
library("IMIFA")
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU <- list()
for (t in 1:nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[t]] <-
    Procrustes(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$U[[t]],
               RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$U[[60001]], translate = TRUE ,dilate = FALSE)$X.new
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}

#------------------------------------------------------------------------------------
## Check complete likelihood for each iteration
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like <- c()
for (t in 1:nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like <-
    c(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like,
      UnDirected_ZIPLPCM_MFM_CompleteLikelihood(X=RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$X[[t]],
                                            U=RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[t]],
                                            beta=RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$beta[t],
                                            nu=RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$nu[[t]],
                                            P=RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP[[t]],
                                            z=RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,],
                                            alpha1=1,alpha2=0.103,omega=0.01,  alpha=3))
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
plot(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like,type = "l",xlab = "",ylab = "", main = "Likelihood",cex.axis = 0.8)

# Check the trace plot of K
plot(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$K,type = "l",xlab = "",ylab = "", main = "K Trace Plot",cex.axis = 0.8)

#------------------------------------------------------------------------------------
# Summarize posterior clustering z by the greedy algorithm proposed by Rastelli and Friel (2018)
# We start from obtaining the marginal posterior mode of the posterior z chain
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States <- c() # initialize a list which will store different clustering states; the clustering states are labeled from 1,2,3... and are put at the 1st,2nd,3rd... row of the matrix, respectively
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration <- list() # store all the t's (iteration number) which provides the same cluster as the clustering state 1,2,3...
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop <- iteration_after_burn_in # the "IterationLoop" which stores all the iteration t's which we focus on, that is, all the iteration t's after burn-in
StatesLabelIndicator = 0 # initialize the label for the clustering states
while (length(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop)!=0){ # if the "IterationLoop" is not empty
  StatesLabelIndicator <- StatesLabelIndicator + 1 # assign the next label to the next clustering state
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState <- RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop[1],] # extract the first clustering state for the "IterationLoop"
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States <- rbind(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States,
                                                                                 RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState) # store the first state within the "IterationLoop" with label "StatesLabelIndicator" in the list which will contain all different unique states
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp <- c() # create a vector to temporarily store all the iteration t's whose clustering is the same as the first clustering state within the "IterationLoop"
  for (t in RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop){ # loop over all the current existing iterations in "IterationLoop"
    if (sum(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,]==RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState)==length(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState)){ # if the t's clustering is the same as the "FirstState"
      RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp <- c(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp,t) # store the iteration t in the temporary vector
    }
  }
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration[[StatesLabelIndicator]] <- RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp # store all the t's as the list element
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop <-
    (iteration_after_burn_in)[-(unlist(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration)-burn_in)] # remove all the iterations we have stored and then move to the next clustering state
}
rownames(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States) <- NULL
nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States) # check the number of different clustering states
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency <- c() # check the number of times one clustering state occurs
for (t in 1:nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States)){
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency <-
    c(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency,
      length(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration[[t]]))
}
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency
sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency,decreasing=TRUE) # check the clustering state frequency in decreasing order
order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency,decreasing=TRUE) # check the corresponding positions
which.max(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency) # find the marginal posterior mode, i.e. the most frequent clustering state
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency[which.max(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency)] # Find the most frequency
library(GreedyEPL) # obtain the summarized z by the greedy algorithm proposed by Rastelli and Friel (2018)
output <- MinimiseEPL(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[iteration_after_burn_in,],
                      list(Kup = 20, loss_type = "VI",decision_init = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States[
                        which.max(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency),]))
output$EPL # check minEVI posterior loss: 0.5421714
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z <- output$decision # Save output as the summarized hat_z

#-----------------------------------------------------------------------------------------
# Obtain summarized U
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_hat_zIteration <- # extract all the iterations whose clustering is identical to hat_z
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration[[
    which(apply(t(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States)==
                  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z,2,sum)==length(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))
  ]]
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_MaxLikeHatz_iteration <- RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_hat_zIteration[
  which.max(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_hat_zIteration])]
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U <-
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_MaxLikeHatz_iteration]]
```

Up to this point, we have obtained the $\hat{\boldsymbol{z}}$ and $\hat{\boldsymbol{U}}$, the interactive 3-dimensional plot of which can be directly produced following:

``` r
library("igraph")
library("RColorBrewer")
library("ggplot2")

My_colors <- c(brewer.pal(10,"RdBu")[c(4,8)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

g_obs <- graph_from_adjacency_matrix(Windsurfers_adj,mode = "undirected",weighted = TRUE)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(Windsurfers_adj))[E(g_obs)$weight]
betw <- betweenness(g_obs)
VertexSize <- sqrt(betw/1.5+mean(betw))*1.5

library("plotly")
fig <- plot_ly() %>%
  add_markers(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(Windsurfers_adj)),
              size=VertexSize,sizes=c(200,400),
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),colors=My_colors[c(1,6,2)]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.3*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "hat_U and hat_z",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

However, it seems to be not easy to find two nice angle to take screenshots of the above interactive 3-d plot in order to show in the **ZIP-LPCM** paper, because the three inferred groups in the above plot lie along the vertical axis.
Considering the fact that our model is invariant under any rotation or translation applied on the latent positions, we propose to apply Procrutes transformation on the $\hat{\boldsymbol{U}}$ to make the three inferred groups lie down along the horizontal plane, i.e., move the latent positions from "portrait" to "landscape":

``` r
library("IMIFA") # translate the latent positions from portrait to landscape
U_ref <- RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U
U_ref[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z==1,] <- c(-3,3,0)
U_ref[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z==2,] <- c(-3,3,0)
U_ref[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z==3,] <- c(3,-3,0)
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U <-
  Procrustes(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U,
             U_ref, translate = TRUE ,dilate = FALSE)$X.new
```

This brings the 3-d plots of the latent positions we illustrate in **Figure 9** of the paper:

``` r
fig <- plot_ly() %>%
  add_markers(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(Windsurfers_adj)),
              size=VertexSize,sizes=c(200,400),
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),colors=My_colors[c(1,6,2)]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.3*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "hat_U and hat_z",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

Such an interactive 3-d plot is uploaded as [`/Interactive 3-d latent positions plots/RDA_Windsurfers_InteractivePlot.html`] of this repository.
The **Figure 9** of the **ZIP-LPCM** paper can thus be reporduced following:

``` r
# Plot the front angle of the latent positions
fig1 <- plot_ly(scene ="scene1") %>% # plot the summarized clustering and hat_U
  add_markers(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(Windsurfers_adj)),
              size=VertexSize,sizes=c(200,400),showlegend = FALSE,
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),colors=My_colors[c(1,6,2)]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig1 <- fig1 %>%
    add_trace(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),showlegend = FALSE,
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.3*E(g_obs)$weight[i]))
}

# Plot the right angle of the latent positions
fig2 <- plot_ly(scene ="scene2") %>% # plot the summarized clustering and hat_U
  add_markers(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(Windsurfers_adj)),
              size=VertexSize,sizes=c(200,400),showlegend = FALSE,
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),colors=My_colors[c(1,6,2)]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig2 <- fig2 %>%
    add_trace(x = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.3*E(g_obs)$weight[i]))
}
Fig1 <- subplot(fig1, fig2)
Fig1 <- Fig1 %>% layout(title = "", margin = list(l = 0,r = 0,b = 0,t = 0,pad = 0),
                        scene = list(domain=list(x=c(0,1/2),y=c(0,1)),
                                     xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                     camera = list(eye = list(x = 1.0, y = 0.6, z = 1.4)),
                                     aspectmode='auto'),
                        scene2 = list(domain=list(x=c(1/2,0.999),y=c(0,1)),
                                      xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                      camera = list(eye = list(x = -0.6, y = 1.0, z = 1.4)),
                                      aspectmode='auto'))
# Fig1
orca(Fig1, "RDA_Windsurfers_hat_U_Y.pdf",scale=1,width=1500,height=700)
```

The rest summary statistics can be obtained by:

``` r
## Summarize beta
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_beta <-
  mean(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$beta[iteration_after_burn_in])
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_beta #  3.614686
#---------------------------------------------------------------------------------------------------------------------------
## Obtain posterior mean of nu, i.e. approximated P_m0
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu <-
  Reduce("+",RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$nu[iteration_after_burn_in])/length(iteration_after_burn_in)
#---------------------------------------------------------------------------------------------------------------------------

library("RColorBrewer")
library("pheatmap")
library("ggplot2")

My_colors <- c(brewer.pal(10,"RdBu")[c(4,8)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

Windsurfers_adj_dataframe <- as.data.frame(Windsurfers_adj)
rownames(Windsurfers_adj_dataframe) <- colnames(Windsurfers_adj_dataframe) <- 1:nrow(Windsurfers_adj)

RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe <- as.data.frame(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu)
rownames(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe) <- colnames(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe) <- 1:nrow(Windsurfers_adj)
#---------------------------------------------------------------------------------------------------------------------------
# Obtain the individual-level p for each iteration
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_p <- list()
library(Rfast)
for (t in 1:nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,],length(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,]),
                  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))==(1:RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))*1
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_p[[t]] <- Z%*%RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP[[t]]%*%t(Z)
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p <- Reduce("+",RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_p[iteration_after_burn_in])/length(iteration_after_burn_in)
diag(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p) <- 0
#---------------------------------------------------------------------------------------------------------------------------
# Obtain the individual-level lambda for each iteration
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_lambda <- list()
library(Rfast)
for (t in 1:nrow(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,],length(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,]),RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))==(1:RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))*1
  Dist_U <- Rfast::Dist(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[t]])
  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_lambda[[t]] <- exp(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1$beta[t]-Dist_U)
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_lambda <- Reduce("+",RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_lambda[iteration_after_burn_in])/length(iteration_after_burn_in)
diag(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_lambda) <- 0
```

And thus the **Figure 10** in the paper can be recovered by:

``` r
annotation_row_z_ref <- as.data.frame(as.factor(matrix(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z,length(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),1)))
colnames(annotation_row_z_ref) <- "z_ref"
annotation_colors_z_ref <- My_colors[c(1,6,2)]
names(annotation_colors_z_ref) <- sort(unique(annotation_row_z_ref$z_ref))
annotation_colors_blank <- rep(brewer.pal(9,"Greys")[1],3)
names(annotation_colors_blank) <- sort(unique(annotation_row_z_ref$z_ref))

RDA_Windsurfers_UnDirected_Y <-
  pheatmap(Windsurfers_adj_dataframe,
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(max(Windsurfers_adj))),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_blank),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE)

RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_Y_hat_z_heatmap <-
  pheatmap(Windsurfers_adj_dataframe[order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(max(Windsurfers_adj))),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))

RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_Non0_heatmap <-
  pheatmap((RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe*(1-dpois(0,RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_lambda)))[order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:8])(10000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))

library("grid")
library("gridExtra")
g <- grid.arrange(RDA_Windsurfers_UnDirected_Y[[4]],
                  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_Y_hat_z_heatmap[[4]],
                  RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_Non0_heatmap[[4]],
                  nrow=1,ncol=3,vp=viewport(width=1, height=1))
Fig <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(colour = "white", linewidth = 0))
Fig
```

And further the heatmap plots of the $\hat{\boldsymbol{p}}$ and $\hat{\boldsymbol{\nu}}$ can be produced following:

``` r
# Plot the hat_nu
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap <-
  pheatmap(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe[order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(10000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap_draw <- cowplot::ggdraw(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap[[4]])+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[1]))
print(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap_draw)

# Plot the hat_p
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe <- as.data.frame(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p)
rownames(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe) <- colnames(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe) <- 1:nrow(Windsurfers_adj)
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap <-
  pheatmap(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe[order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:3])(10000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap_draw <- cowplot::ggdraw(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap[[4]])+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[1]))
print(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap_draw)
```

### 2.3 Train Bombing Network

Similar to the **Windsurfers** real network, the **Train bombing** network is also an undirected network and reference clustering or exogenous node attributes are not available to us.
Thus an **unsupervised ZIP-LPCM** implementation is applied in practice for this network along with the post processing steps:

``` r
# Train Bombing undirected real network unsupervised ZIP-LPCM T = 60000 round 1
set.seed(1)
start.time <- Sys.time()
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1 <- 
  MwG_UnDirected_ZIPLPCM(Y = Train_bombing_adj,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=9,
                         sigma2prop_beta=0.3^2,sigma2prop_U=0.7,d=3,z=1:nrow(Train_bombing_adj),
                         p_eject=0.5)
end.time <- Sys.time()
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_time <- end.time - start.time
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_time
# Time difference of 1.747663 hours

#--------------------------------------------------------------------------------------------------------------------------
# Post-process the output

# Define the burn in
iteration_after_burn_in <- 30002:60001
burn_in <- 30001

## check U acceptance rate
apply(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)
mean(apply(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)) # 0.3209318
## check beta acceptance rate
mean(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$acceptance_count_beta[iteration_after_burn_in-1]) # 0.2292

# Apply label switching on the post clustering z
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz <- matrix(NA,nrow=nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z),ncol=ncol(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z))
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP <- list()
for (t in 1:nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  LS_temp <- LabelSwitching_SG2003(z = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z[t,],
                                   matrix_KbyK = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$P[[t]])
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,] <- LS_temp$z
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP[[t]] <- LS_temp$matrix_KbyK
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
# Apply Procrustes Transform on posterior latent positions U
library("IMIFA")
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU <- list()
for (t in 1:nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[t]] <-
    Procrustes(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$U[[t]],
               RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$U[[60001]], translate = TRUE ,dilate = FALSE)$X.new
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}

#--------------------------------------------------------------------------------------------------------------------------
## Check complete likelihood for each iteration
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like <- c()
for (t in 1:nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like <-
    c(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like,
      UnDirected_ZIPLPCM_MFM_CompleteLikelihood(X=RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$X[[t]],
                                            U=RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[t]],
                                            beta=RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$beta[t],
                                            nu=RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$nu[[t]],
                                            P=RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP[[t]],
                                            z=RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,],
                                            alpha1=1,alpha2=0.103,omega=0.01,  alpha=3))
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
plot(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like,type = "l",xlab = "",ylab = "", main = "Likelihood",cex.axis = 0.8)

# Check the trace plot of K
plot(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$K,type = "l",xlab = "",ylab = "", main = "K Trace Plot",cex.axis = 0.8)

#--------------------------------------------------------------------------------------------------------------------------
# Summarize posterior clustering z by the greedy algorithm proposed by Rastelli and Friel (2018)
# We start from obtaining the marginal posterior mode of the posterior z chain
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States <- c() # initialize a list which will store different clustering states; the clustering states are labeled from 1,2,3... and are put at the 1st,2nd,3rd... row of the matrix, respectively
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration <- list() # store all the t's (iteration number) which provides the same cluster as the clustering state 1,2,3...
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop <- iteration_after_burn_in # the "IterationLoop" which stores all the iteration t's which we focus on, that is, all the iteration t's after burn-in
StatesLabelIndicator = 0 # initialize the label for the clustering states
while (length(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop)!=0){ # if the "IterationLoop" is not empty
  StatesLabelIndicator <- StatesLabelIndicator + 1 # assign the next label to the next clustering state
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState <- RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop[1],] # extract the first clustering state for the "IterationLoop"
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States <- rbind(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States,
                                                                                 RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState) # store the first state within the "IterationLoop" with label "StatesLabelIndicator" in the list which will contain all different unique states
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp <- c() # create a vector to temporarily store all the iteration t's whose clustering is the same as the first clustering state within the "IterationLoop"
  for (t in RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop){ # loop over all the current existing iterations in "IterationLoop"
    if (sum(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,]==RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState)==length(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_FirstState)){ # if the t's clustering is the same as the "FirstState"
      RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp <- c(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp,t) # store the iteration t in the temporary vector
    }
  }
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration[[StatesLabelIndicator]] <- RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration_temp # store all the t's as the list element
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_IterationLoop <-
    (iteration_after_burn_in)[-(unlist(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration)-burn_in)] # remove all the iterations we have stored and then move to the next clustering state
}
rownames(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States) <- NULL
nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States) # check the number of different clustering states
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency <- c() # check the number of times one clustering state occurs
for (t in 1:nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States)){
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency <-
    c(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency,
      length(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration[[t]]))
}
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency
sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency,decreasing=TRUE) # check the clustering state frequency in decreasing order
order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency,decreasing=TRUE) # check the corresponding positions
which.max(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency) # find the marginal posterior mode, i.e. the most frequent clustering state
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency[which.max(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency)] # Find the most frequency
library(GreedyEPL) # obtain the summarized z by the greedy algorithm proposed by Rastelli and Friel (2018)
output <- MinimiseEPL(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[iteration_after_burn_in,],
                      list(Kup = 20, loss_type = "VI",decision_init = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States[
                        which.max(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesFrequency),]))
output$EPL # check minEVI posterior loss: 0.4427744
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z <- output$decision # Save output as summarized hat_z

#--------------------------------------------------------------------------------------------------------------------------
# Obtain summarized U
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_hat_zIteration <- # extract all the iterations whose clustering is identical to hat_z
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_StatesIteration[[
    which(apply(t(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_States)==
                  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z,2,sum)==length(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))
  ]]
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_MaxLikeHatz_iteration <- RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_hat_zIteration[
  which.max(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_Like[RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz_hat_zIteration])]
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U <-
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_MaxLikeHatz_iteration]]
```

After obtaining the summarized $\hat{\boldsymbol{z}}$ and $\hat{\boldsymbol{U}}$, the corresponding 3-d interactive plot can be produced following:

``` r
library("igraph")
library("RColorBrewer")
library("ggplot2")

My_colors <- c(brewer.pal(10,"RdBu")[c(4,8)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

g_obs <- graph_from_adjacency_matrix(Train_bombing_adj,mode = "undirected",weighted = TRUE)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(Train_bombing_adj))[E(g_obs)$weight]
betw <- betweenness(g_obs)
VertexSize <- sqrt(betw/1.5+mean(betw))*1

library("plotly")
fig <- plot_ly() %>%
  add_markers(x = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,1],
              y = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,2],
              z = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(Train_bombing_adj)),
              size=VertexSize,sizes=c(200,400),
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),colors=My_colors[c(6,2,8,9)]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.75*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "hat_U and hat_z",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

Such a 3-d plot can be downloaded at [`/Interactive 3-d latent positions plots/RDA_TrainBombing_InteractivePlot.html`] of this repository.
Following the above code for the 3-d plotting, the **Figure 11** in the **ZIP-LPCM** paper can be reproduced by:

``` r
# Plot the front angle of the latent positions
fig1 <- plot_ly(scene ="scene1") %>%
  add_markers(x = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,1],
              y = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,2],
              z = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(Train_bombing_adj)),
              size=VertexSize,sizes=c(200,400),showlegend = FALSE,
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),colors=My_colors[c(6,2,8,9)]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig1 <- fig1 %>%
    add_trace(x = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.75*E(g_obs)$weight[i]))
}


# Plot the right angle of the latent positions
fig2 <- plot_ly(scene ="scene2") %>%
  add_markers(x = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,1],
              y = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,2],
              z = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(Train_bombing_adj)),
              size=VertexSize,sizes=c(200,400),showlegend = FALSE,
              stroke = I("black"), span = I(1),
              color=as.factor(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),colors=My_colors[c(6,2,8,9)]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig2 <- fig2 %>%
    add_trace(x = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.75*E(g_obs)$weight[i]))
}

Fig1 <- subplot(fig1, fig2)
Fig1 <- Fig1 %>% layout(title = "", margin = list(l = 0,r = 0,b = 0,t = 0,pad = 0),
                        scene = list(domain=list(x=c(0,1/2),y=c(0,1)),
                                     xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                     camera = list(eye = list(x = -1.1, y = -1.5, z = 0.85)),
                                     aspectmode='auto'),
                        scene2 = list(domain=list(x=c(1/2,0.999),y=c(0,1)),
                                      xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                      camera = list(eye = list(x = 1.1, y = -1.5, z = 0.85)),
                                      aspectmode='auto'))
# Fig1
orca(Fig1, "RDA_TrainBombing_hat_U_Y.pdf",scale=1,width=1800,height=850)
```

Finally, the rest summary statistics and the corresponding heatmap plots as well as the **Figure 11** in the paper can be produced following the code below:

``` r
## Summarize beta
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_beta <-
  mean(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$beta[iteration_after_burn_in])
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_beta #  1.057638
#--------------------------------------------------------------------------------------------------------------------------
## Obtain posterior mean of nu, i.e. approximated P_m0
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu <-
  Reduce("+",RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$nu[iteration_after_burn_in])/length(iteration_after_burn_in)
#--------------------------------------------------------------------------------------------------------------------------
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
My_colors <- c(brewer.pal(10,"RdBu")[c(4,8)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

Train_bombing_adj_dataframe <- as.data.frame(Train_bombing_adj)
rownames(Train_bombing_adj_dataframe) <- colnames(Train_bombing_adj_dataframe) <- 1:nrow(Train_bombing_adj)

RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe <- as.data.frame(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu)
rownames(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe) <- colnames(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe) <- 1:nrow(Train_bombing_adj)
#--------------------------------------------------------------------------------------------------------------------------
# Obtain p for each iteration
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_p <- list()
library(Rfast)
for (t in 1:nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,],length(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,]),
                  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))==(1:RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))*1
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_p[[t]] <- Z%*%RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSP[[t]]%*%t(Z)
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p <- Reduce("+",RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_p[iteration_after_burn_in])/length(iteration_after_burn_in)
diag(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p) <- 0

#--------------------------------------------------------------------------------------------------------------------------
# Obtain lambda for each iteration
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_lambda <- list()
library(Rfast)
for (t in 1:nrow(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,],length(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_LSz[t,]),RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))==(1:RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$K[t]))*1
  Dist_U <- Rfast::Dist(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_PTU[[t]])
  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_lambda[[t]] <- exp(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1$beta[t]-Dist_U)
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_lambda <- Reduce("+",RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_lambda[iteration_after_burn_in])/length(iteration_after_burn_in)
diag(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_lambda) <- 0

#--------------------------------------------------------------------------------------------------------------------------
# Make the heatmap plots
annotation_row_z_ref <- as.data.frame(as.factor(matrix(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z,length(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),1)))
colnames(annotation_row_z_ref) <- "z_ref"
annotation_colors_z_ref <- My_colors[c(6,2,8,9)]
names(annotation_colors_z_ref) <- sort(unique(annotation_row_z_ref$z_ref))
annotation_colors_blank <- rep(brewer.pal(9,"Greys")[1],4)
names(annotation_colors_blank) <- sort(unique(annotation_row_z_ref$z_ref))

RDA_TrainBombing_UnDirected_Y <-
  pheatmap(Train_bombing_adj_dataframe,
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(max(Train_bombing_adj))),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_blank),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE)

RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_Y_hat_z_heatmap <-
  pheatmap(Train_bombing_adj_dataframe[order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(max(Train_bombing_adj))),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))


RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_Non0_heatmap <-
  pheatmap((RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe*(1-dpois(0,RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_lambda)))[order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:5])(2000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))

library("grid")
library("gridExtra")
g <- grid.arrange(RDA_TrainBombing_UnDirected_Y[[4]],
                  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_Y_hat_z_heatmap[[4]],
                  RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_Non0_heatmap[[4]],
                  nrow=1,ncol=3,vp=viewport(width=1, height=1))
Fig <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(colour = "white", linewidth = 0))
Fig

#--------------------------------------------------------------------------------------------------------------------------
# Plot the hat_nu
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap <-
  pheatmap(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_nu_dataframe[order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:6])(2000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap_draw <- cowplot::ggdraw(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap[[4]])+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[1]))
print(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_P_m0_heatmap_draw)
# Plot the hat_p
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe <- as.data.frame(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p)
rownames(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe) <- colnames(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe) <- 1:nrow(Train_bombing_adj)
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap <-
  pheatmap(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_dataframe[order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z),order(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:6])(100)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_z))!=0)))
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap_draw <- cowplot::ggdraw(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap[[4]])+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[1]))
print(RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_T60k_R1_hat_p_heatmap_draw)
```

### 2.4 Summit Co-attendance Criminality Network

The last real network we work on is the **Summit Co-attendance Criminality** network.
This network is an undirected network and exogenous node attributes are available to us, i.e., `RDA_criminalNet$A` which is available in Section 1.3 of this tutorial.
The implementations and the post-processing follows:

``` r
# Summit co-attendance criminality undirected real network supervised ZIP-LPCM T = 60000 round 1
set.seed(1)
start.time <- Sys.time()
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1 <- 
  MwG_UnDirected_ZIPLPCM(Y = RDA_criminalNet$Y,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=6,
                         sigma2prop_beta=0.125^2,sigma2prop_U=0.5,d=3,z=1:nrow(RDA_criminalNet$Y),
                         p_eject=0.5,A=RDA_criminalNet$A,omega_c=1)
end.time <- Sys.time()
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_time <- end.time - start.time
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_time
# Time difference of 2.744616 hours

#--------------------------------------------------------------------------------------------------------------------------
# Post-process the output

# Define the burn in
iteration_after_burn_in <- 30002:60001
burn_in <- 30001

## check U acceptance rate
apply(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)
mean(apply(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$acceptance_count_U[iteration_after_burn_in-1,],2,mean)) # 0.2191056
## check beta acceptance rate
mean(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$acceptance_count_beta[iteration_after_burn_in-1]) # 0.2617

# Apply label switching on the post clustering z
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz <- matrix(NA,nrow=nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z),ncol=ncol(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z))
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSP <- list()
for (t in 1:nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z)){
  LS_temp <- LabelSwitching_SG2003(z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z[t,],
                                   matrix_KbyK = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$P[[t]])
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[t,] <- LS_temp$z
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSP[[t]] <- LS_temp$matrix_KbyK
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
# Apply Procrustes Transform on posterior latent positions U
library("IMIFA")
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_PTU <- list()
for (t in 1:nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z)){
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_PTU[[t]] <-
    Procrustes(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$U[[t]],
               RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$U[[60001]], translate = TRUE ,dilate = FALSE)$X.new
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}

#--------------------------------------------------------------------------------------------------------------------------
## Check complete likelihood for each iteration
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_Like <- c()
for (t in 1:nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z)){
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_Like <-
    c(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_Like,
      UnDirected_ZIPLPCM_MFM_CompleteLikelihood(X=RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$X[[t]],
                                            U=RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_PTU[[t]],
                                            beta=RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$beta[t],
                                            nu=RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$nu[[t]],
                                            P=RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSP[[t]],
                                            z=RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[t,],
                                            alpha1=1,alpha2=0.103,omega=0.01,  alpha=3,
                                            A=RDA_criminalNet$A,omega_c=1))
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
plot(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_Like,type = "l",xlab = "",ylab = "", main = "Likelihood",cex.axis = 0.8)

# Check the trace plot of K
plot(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$K,type = "l",xlab = "",ylab = "", main = "K Trace Plot",cex.axis = 0.8)

#--------------------------------------------------------------------------------------------------------------------------
# Summarize posterior clustering z by the greedy algorithm proposed by Rastelli and Friel (2018)
# We start from obtaining the marginal posterior mode of the posterior z chain
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States <- c() # initialize a list which will store different clustering states; the clustering states are labeled from 1,2,3... and are put at the 1st,2nd,3rd... row of the matrix, respectively
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration <- list() # store all the t's (iteration number) which provides the same cluster as the clustering state 1,2,3...
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop <- iteration_after_burn_in # the "IterationLoop" which stores all the iteration t's which we focus on, that is, all the iteration t's after burn-in
StatesLabelIndicator = 0 # initialize the label for the clustering states
while (length(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop)!=0){ # if the "IterationLoop" is not empty
  StatesLabelIndicator <- StatesLabelIndicator + 1 # assign the next label to the next clustering state
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_FirstState <- RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop[1],] # extract the first clustering state for the "IterationLoop"
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States <- rbind(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States,
                                                                              RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_FirstState) # store the first state within the "IterationLoop" with label "StatesLabelIndicator" in the list which will contain all different unique states
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp <- c() # create a vector to temporarily store all the iteration t's whose clustering is the same as the first clustering state within the "IterationLoop"
  for (t in RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop){ # loop over all the current existing iterations in "IterationLoop"
    if (sum(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[t,]==RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_FirstState)==length(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_FirstState)){ # if the t's clustering is the same as the "FirstState"
      RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp <- c(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp,t) # store the iteration t in the temporary vector
    }
  }
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration[[StatesLabelIndicator]] <- RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration_temp # store all the t's as the list element
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_IterationLoop <-
    (iteration_after_burn_in)[-(unlist(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration)-burn_in)] # remove all the iterations we have stored and then move to the next clustering state
}
rownames(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States) <- NULL
nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States) # check the number of different clustering states
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency <- c() # check the number of times one clustering state occurs
for (t in 1:nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States)){
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency <-
    c(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency,
      length(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration[[t]]))
}
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency
sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency,decreasing=TRUE) # check the clustering state frequency in decreasing order
order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency,decreasing=TRUE) # check the corresponding positions
which.max(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency) # find the marginal posterior mode, i.e. the most frequent clustering state
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency[which.max(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency)] # Find the most frequency
library(GreedyEPL) # obtain the summarized z by the greedy algorithm proposed by Rastelli and Friel (2018)
output <- MinimiseEPL(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,],
                      list(Kup = 20, loss_type = "VI",decision_init = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States[
                        which.max(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesFrequency),]))
table(output$decision,RDA_criminalNet$OurA,dnn = c("","")) # check whether the the output is the same as the marginal posterior mode
output$EPL # check minEVI posterior loss: 0.2661251
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z <- output$decision # Save output as summarized z


#--------------------------------------------------------------------------------------------------------------------------
# Obtain summarized U
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_hat_zIteration <- # extract all the iterations whose clustering is identical to hat_z
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_StatesIteration[[
    which(apply(t(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_States)==
                  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z,2,sum)==length(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z))
  ]]
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_MaxLikeHatz_iteration <- RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_hat_zIteration[
  which.max(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_Like[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz_hat_zIteration])]
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U <-
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_PTU[[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_MaxLikeHatz_iteration]]
```

Provided with the reference clustering $`\boldsymbol{z}^*`$, i.e., the role-affiliation information stored as `RDA_criminalNet$RoleAffiliation`, as well as the summarized latent positions $\hat{\boldsymbol{U}}$, we first obtain the 3-dimensional interactive plot of $\hat{\boldsymbol{U}}$ plotting along with the $`\boldsymbol{z}^*`$.
This is the one we illustrate as the 1st row plots of **Figure 13** in the paper:

``` r
library("igraph")
library("RColorBrewer")
library("ggplot2")

My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

g_obs <- graph_from_adjacency_matrix(RDA_criminalNet$Y,mode = "undirected",weighted = TRUE)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(RDA_criminalNet$Y))[E(g_obs)$weight]
betw <- betweenness(g_obs)
VertexSize <- sqrt(betw/1.5+mean(betw))*1.25

library("plotly")
# Plot the reference clustering z* and hat_U
fig <- plot_ly() %>%
  add_markers(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(RDA_criminalNet$Y),"<br>z*:",RDA_criminalNet$RoleAffiliation,"<br>c:",RDA_criminalNet$A),
              size=VertexSize,sizes=c(200,400),
              color=as.factor(RDA_criminalNet$RoleAffiliation),colors=My_colors[1:10],
              stroke = I("black"), span = I(1),
              symbol=as.factor(RDA_criminalNet$role),symbols = c("circle","square")
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.5*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "hat_U and z*",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

The 3-d interactive plot of the $\hat{\boldsymbol{U}}$ and $`\boldsymbol{z}^*`$ is available at [`/Interactive 3-d latent positions plots/RDA_CriminalSummit_InteractivePlot_Ref_z.html`] of this repository.
Different from the 3-d interactive plots we show in the previous sections, if the readers put the mouse pointer on each node of the interactive plot, the comment bracket will also show the exogenous node attributes $\boldsymbol{c}$ we used in practice in our experiments.
Note that the **Windsurfers**, **Train Bombing** and the **Summit Co-attendance Criminality** networks are undirected networks, the comment bracket of each edge in the 3-d interactive plots no longer show the direction of the edge.

Similarly, the 3-d interactive plot of the $\hat{\boldsymbol{U}}$ and $\hat{\boldsymbol{z}}$ is available to be downloaded at [`/Interactive 3-d latent positions plots/RDA_CriminalSummit_InteractivePlot_Hat_z.html`] of this repository, and can be reproduced following the code:

``` r
# Plot the hat_z and hat_U
fig <- plot_ly() %>%
  add_markers(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(RDA_criminalNet$Y),"<br>z*:",RDA_criminalNet$RoleAffiliation,"<br>c:",RDA_criminalNet$A),
              size=VertexSize,sizes=c(200,400),
              color=as.factor(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z),colors=My_colors[c(3,1,6,9,4,5,7,10,8,2)],
              stroke = I("black"), span = I(1),
              symbol=as.factor(RDA_criminalNet$role),symbols = c("circle","square")
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.5*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "hat_U and hat_z",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

As we discussed in the **ZIP-LPCM** paper, it's interesting that there are three heterogeneous core criminal suspects (bosses) whose corresponding inferred latent positions are close to each other and to the orange affiliation.
These three bosses are, respectively, a blue boss node no.53, an orange boss node no.41 and a green boss node no.38.
Thus we can check within how many posterior samples of clustering each pair of these three bosses are clustered together in one group:

``` r
### Check clustering of the heterogeneous bosses after burn-in
## Check how often orange boss 41 and blue boss 53 in the same group
sum((RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,41]==
       RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,53]))/length(iteration_after_burn_in)
# 0.6912667
## Check how often green boss 38 and blue boss 53 in the same group
sum((RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,38]==
       RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,53]))/length(iteration_after_burn_in)
# 0.2982667
## Check how often green boss 38 and orange boss 41 in the same group
sum((RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,38]==
       RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,41]))/length(iteration_after_burn_in)
# 0.002133333
## Check how often green boss 38 and orange boss 41 and blue boss 53 in the same group
sum((RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,38]==
       RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,41])*
      (RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,38]==
         RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,53]))/length(iteration_after_burn_in)
# 0.002133333
```

Thus it shows that the three heterogeneous bosses are not likely to be clustered into the same group.
Though the clustering of the blue boss node no.53 is arguable as shown above, this boss is more likely to be clustered with the orange boss node no.41 who is clustered inside the core group of the orange affiliation based on our inferred clustering $\hat{\boldsymbol{z}}$.
The green boss node no.38 is not likely to be within the core of the orange affiliation and most of the time is clustered with another orange boss node no.35 in the same group:

``` r
## Check how often another orange boss 35 and green boss 38 in the same group
sum((RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,38]==
       RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[iteration_after_burn_in,35]))/length(iteration_after_burn_in)
# 0.9975
```

The **Figure 13** illustrated in the **ZIP-LPCM** paper can be reproduced by:

``` r
# First plot the reference clustering z* and hat_U
# Plot the front angle of the latent positions
fig1 <- plot_ly(scene ="scene1") %>%
  add_markers(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(RDA_criminalNet$Y),"<br>z*:",RDA_criminalNet$RoleAffiliation,"<br>c:",RDA_criminalNet$A),
              size=VertexSize,sizes=c(200,400), showlegend = FALSE,
              color=as.factor(RDA_criminalNet$RoleAffiliation),colors=My_colors[1:10],
              stroke = I("black"), span = I(1),
              symbol=as.factor(RDA_criminalNet$role),symbols = c("circle","square")
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig1 <- fig1 %>%
    add_trace(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.5*E(g_obs)$weight[i]))
}

# Plot the right angle of the latent positions
fig2 <- plot_ly(scene ="scene2") %>%
  add_markers(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(RDA_criminalNet$Y),"<br>z*:",RDA_criminalNet$RoleAffiliation,"<br>c:",RDA_criminalNet$A),
              size=VertexSize,sizes=c(200,400),showlegend = FALSE,
              color=as.factor(RDA_criminalNet$RoleAffiliation),colors=My_colors[1:10],
              stroke = I("black"), span = I(1),
              symbol=as.factor(RDA_criminalNet$role),symbols = c("circle","square")
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig2 <- fig2 %>%
    add_trace(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.5*E(g_obs)$weight[i]))
}

Fig1 <- subplot(fig1, fig2) 
Fig1 <- Fig1 %>% layout(title = "", margin = list(l = 0,r = 0,b = 0,t = 0,pad = 0),
                        scene = list(domain=list(x=c(0,1/2),y=c(0,1)),
                                     xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                     camera = list(eye = list(x = 0, y = 1.50, z = 0.50)),
                                     aspectmode='auto'),
                        scene2 = list(domain=list(x=c(1/2,0.999),y=c(0,1)),
                                      xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                      camera = list(eye = list(x = -1.50*sqrt(3)/2, y = 1.50/2, z = 0.50)),
                                      aspectmode='auto'))
# Fig1
orca(Fig1, "RDA_CriminalSummit_hat_U_RoleAffiliation.pdf",scale=1,width=1800,height=850)

#--------------------------------------------------------------------------------------------------------------------------
# Then plot the summarized clustering hat_z and hat_U

# Plot the front angle of the latent positions
fig3 <- plot_ly(scene ="scene3") %>%
  add_markers(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(RDA_criminalNet$Y),"<br>z*:",RDA_criminalNet$RoleAffiliation,"<br>c:",RDA_criminalNet$A),
              size=VertexSize,sizes=c(200,400),showlegend = FALSE,
              color=as.factor(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z),colors=My_colors[c(3,1,6,9,4,5,7,10,8,2)],
              stroke = I("black"), span = I(1),
              symbol=as.factor(RDA_criminalNet$role),symbols = c("circle","square")
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig3 <- fig3 %>%
    add_trace(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.5*E(g_obs)$weight[i]))
}

# Plot the right angle of the latent positions
fig4 <- plot_ly(scene ="scene4") %>% # plot the hat_z and hat_U
  add_markers(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[,3],
              text=paste("Node:",1:nrow(RDA_criminalNet$Y),"<br>z*:",RDA_criminalNet$RoleAffiliation,"<br>c:",RDA_criminalNet$A),
              size=VertexSize,sizes=c(200,400),showlegend = FALSE,
              color=as.factor(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z),colors=My_colors[c(3,1,6,9,4,5,7,10,8,2)],
              stroke = I("black"), span = I(1),
              symbol=as.factor(RDA_criminalNet$role),symbols = c("circle","square")
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig4 <- fig4 %>%
    add_trace(x = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],1],
              y = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],2],
              z = RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.5*E(g_obs)$weight[i]))
}
Fig2 <- subplot(fig3, fig4) 
Fig2 <- Fig2 %>% layout(title = "", margin = list(l = 0,r = 0,b = 0,t = 0,pad = 0),
                        scene3 = list(domain=list(x=c(0,1/2),y=c(0,1)),
                                      xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                      camera = list(eye = list(x = 0, y = 1.50, z = 0.50)),
                                      aspectmode='auto'),
                        scene4 = list(domain=list(x=c(1/2,0.999),y=c(0,1)),
                                      xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                      camera = list(eye = list(x = -1.50*sqrt(3)/2, y = 1.50/2, z = 0.50)),
                                      aspectmode='auto'))
# Fig2
orca(Fig2, "RDA_CriminalSummit_hat_U_hat_z.pdf",scale=1,width=1800,height=850)
```

Similar to all the previous experiments, we can follow the code below to obtain the summary statistics: (i) the posterior mean of intercept $\hat{\beta}$; (ii) the posterior mean of unusual zero indicator $\hat{\boldsymbol{\nu}}$ which approximates the conditional probability of missing zeros provided that the corresponding interactions are zeros, i.e., the `P_m0` we denote in the code and the equation (22) of the **ZIP-LPCM** paper; (iii) the individual-level unusual probability $\boldsymbol{p}$ which is a $N \times N$ matrix with each entry $i,j$ being the $p_{z_iz_j}$ for each posterior state, and thus the posterior mean brings the $\hat{\boldsymbol{p}}$ accounting for the uncertainty of the posterior clustering; (iv) the individual-level Poisson rate $N \times N$ matrix `lambda` we denote in the code with each entry being the $\lambda_{ij}=\text{exp}(\beta-||\boldsymbol{u_i}-\boldsymbol{u_j}||)$ obtained by the corresponding $\beta$ and $\boldsymbol{U}$ for each posterior state, and the corresponding posterior mean `hat_lambda` can thus be obtained.

``` r
## Summarize beta
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_beta <-
  mean(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$beta[iteration_after_burn_in])
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_beta # 2.169797
#--------------------------------------------------------------------------------------------------------------------------
## Obtain posterior mean of nu, i.e. approximated P_m0
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_nu <-
  Reduce("+",RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$nu[iteration_after_burn_in])/length(iteration_after_burn_in)
#--------------------------------------------------------------------------------------------------------------------------
# Obtain individual-level unusual zero probability p for each iteration
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_p <- list()
library(Rfast)
for (t in 1:nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[t,],length(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[t,]),
                  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$K[t]))==(1:RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$K[t]))*1
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_p[[t]] <- Z%*%RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSP[[t]]%*%t(Z)
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p <- Reduce("+",RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_p[iteration_after_burn_in])/length(iteration_after_burn_in)

#--------------------------------------------------------------------------------------------------------------------------
# Obtain individual-level Poisson rate lambda for each iteration
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_lambda <- list()
library(Rfast) # for Rfast::Dist()
for (t in 1:nrow(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$z)){
  Z <- t(t(matrix(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[t,],length(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_LSz[t,]),RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$K[t]))==(1:RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$K[t]))*1
  Dist_U <- Rfast::Dist(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_PTU[[t]])
  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_lambda[[t]] <- exp(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1$beta[t]-Dist_U)
  if ((t%%1000) == 0){cat("t=",t,"\n")}
}
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_lambda <- Reduce("+",RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_lambda[iteration_after_burn_in])/length(iteration_after_burn_in)
```

The heatmap plots shown in **Figure 14** of the **ZIP-LPCM** paper can be recoverred by first creating a new clustering which label-switches the $\hat{\boldsymbol{z}}$ so that the inferred groups which contain bosses are sorted to have higher group numbers, that is,

``` r
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual <- RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==2]<-4
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==10]<-2
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==1]<-3
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==5]<-1
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==6]<-5
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==3]<-9
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==7]<-7
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==9]<-8
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==4]<-6
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z==8]<-10
```

In this way, the left color bars of **Figure 14** plots showing the reference clustering $\boldsymbol{z}^*$ would let the groups which contain bosses sit in the bottom and let the groups only containing affiliates sit at top.
Similarly for the top color bars.
Note that this rearragement would not affect any results we obtained, but it can help us have a clearer view of how the core-periphery architecture is structured by criminals based on the heatmap plots of the adcacency matrices, inferred conditional unusual zero probability and so on.
Thus the heatmap plots in **Figure 14** can be reproduced via:

``` r
library("RColorBrewer")
library("pheatmap")
library("ggplot2")

My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

RDA_criminalNet_Y_dataframe <- as.data.frame(RDA_criminalNet$Y)
rownames(RDA_criminalNet_Y_dataframe) <- colnames(RDA_criminalNet_Y_dataframe) <- 1:nrow(RDA_criminalNet$Y)

RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe <- as.data.frame(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_nu)
rownames(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe) <- colnames(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe) <- 1:nrow(RDA_criminalNet$Y)

annotation_row_z_ref <- as.data.frame(as.factor(matrix(RDA_criminalNet$RoleAffiliation,length(RDA_criminalNet$RoleAffiliation),1)))
colnames(annotation_row_z_ref) <- "z_ref"
annotation_colors_z_ref <- My_colors[c(1:10)]
names(annotation_colors_z_ref) <- sort(unique(annotation_row_z_ref$z_ref))

RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_Y_hat_z_heatmap <-
  pheatmap(RDA_criminalNet_Y_dataframe[order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual),order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd"))(max(RDA_criminalNet$Y))),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual))!=0)),gaps_col=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual))!=0))
  )

RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_P_m0_heatmap <-
  pheatmap(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe[order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual),order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:8])(5000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual))!=0)),gaps_col=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual))!=0))
  )

RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_P_m0_Non0_heatmap <-
  pheatmap((RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_nu_dataframe*(1-dpois(0,RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_Lambdaij)))[order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual),order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:6])(5000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual))!=0)),gaps_col=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z_Manual))!=0))
  )

library("grid")
library("gridExtra")
g <- grid.arrange(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_Y_hat_z_heatmap[[4]],
                  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_P_m0_heatmap[[4]],
                  RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_P_m0_Non0_heatmap[[4]],
                  nrow=1,ncol=3,vp=viewport(width=1, height=1))
Fig <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(colour = "white", linewidth = 0))
Fig
```

Similarly, we can also produce the heatmap plot of the summarized individual-level unusual zero probability $\hat{\boldsymbol{p}}$ here:

``` r
# Plot the posterior mean summary statistic hat_p
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe <- as.data.frame(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p)
rownames(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe) <- colnames(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe) <- 1:nrow(RDA_criminalNet$Y)
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap <-
  pheatmap(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_dataframe[order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z),order(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z)],
           color=c(brewer.pal(9,"Greys")[3],colorRampPalette(brewer.pal(9,"YlOrRd")[1:7])(5000)),cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)),gaps_col=c(which(diff(sort(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_z))!=0)))
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap_draw <- cowplot::ggdraw(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap[[4]])+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[3]))
print(RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_T60k_R1_hat_p_heatmap_draw)
```

### 2.5 Comparisons to 2-Dimensional Inference

This subsection corresponds to **Section 5.5** of the **ZIP-LPCM** paper where we compare the 3-d inference output to the corresponding 2-d inference output obtained following the similar steps shown in the previous subsections of this tutorial.
The only difference in the 2-d inference is that the dimension of the latent positions is set as $d=2$ during the implementations.
Here we provide some example code to begin with and omit the rest repeated post-processing code.
The readers can simply modify the code variable names of 3-d cases for the implementations of 2-d inference.

``` r
# Sampson monks directed real network supervised ZIP-LPCM d=2 T = 60000 round 1 
set.seed(1)
start.time <- Sys.time()
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_2d_T60k_R1 <- 
  MwG_Directed_ZIPLPCM(Y = SampsonMonks_Directed_adj,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=9,
                       sigma2prop_beta=0.4^2,sigma2prop_U=0.2,d=2,z=1:nrow(SampsonMonks_Directed_adj),
                       p_eject=0.5,A=sampson_monks_group_cloisterville_NodeA,omega_c=1)
end.time <- Sys.time()
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_2d_T60k_R1_time <- end.time - start.time
RDA_SampsonMonks_Directed_ZIPLPCM_Sup_2d_T60k_R1_time
# Time difference of 34.87837 mins


#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# Windsurfers undirected real network unsupervised ZIP-LPCM d=2 T = 60000 round 1
set.seed(1)
start.time <- Sys.time()
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_2d_T60k_R1 <- 
  MwG_UnDirected_ZIPLPCM(Y = Windsurfers_adj,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=19,
                         sigma2prop_beta=0.15^2,sigma2prop_U=0.2,d=2,z=1:nrow(Windsurfers_adj),
                         p_eject=0.5)
end.time <- Sys.time()
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_2d_T60k_R1_time <- end.time - start.time
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_2d_T60k_R1_time
# Time difference of 1.198135 hours

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# Train Bombing undirected real network unsupervised ZIP-LPCM d=2 T = 60000 round 1
set.seed(1)
start.time <- Sys.time()
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_2d_T60k_R1 <- 
  MwG_UnDirected_ZIPLPCM(Y = Train_bombing_adj,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=9,
                         sigma2prop_beta=0.3^2,sigma2prop_U=0.7,d=2,z=1:nrow(Train_bombing_adj),
                         p_eject=0.5)
end.time <- Sys.time()
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_2d_T60k_R1_time <- end.time - start.time
RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_2d_T60k_R1_time
# Time difference of 1.570117 hours

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

# Summit co-attendance criminality undirected real network supervised ZIP-LPCM d=2 T = 60000 round 1
set.seed(1)
start.time <- Sys.time()
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_2d_T60k_R1 <- 
  MwG_UnDirected_ZIPLPCM(Y = RDA_criminalNet$Y,T = 60000,omega=0.01,alpha1=1,alpha2=0.103,alpha=3,beta1=1,beta2=6,
                         sigma2prop_beta=0.125^2,sigma2prop_U=0.5,d=2,z=1:nrow(RDA_criminalNet$Y),
                         p_eject=0.5,A=RDA_criminalNet$A,omega_c=1)
end.time <- Sys.time()
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_2d_T60k_R1_time <- end.time - start.time
RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_2d_T60k_R1_time
# Time difference of 2.573009 hours
```

However, instead of the plots of 3-dimensional latent positions, we need some new code for visualizing the 2-dimensional latent positions for each real network as shown below.

``` r
# 2-dimensional visualization of Sampson monks
library("igraph")
library("RColorBrewer")
My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])
g_obs_S <- graph_from_adjacency_matrix(SampsonMonks_Directed_adj,mode = "directed",weighted = TRUE)
E(g_obs_S)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(SampsonMonks_Directed_adj))[E(g_obs_S)$weight]
betw <- betweenness(g_obs_S)
V(g_obs_S)$size <- sqrt(betw/1.5+mean(betw))*0.75 # set the vertex size
V(g_obs_S)$frame.color <- "black"
V(g_obs_S)$label <- "" 
Customized_colors <- My_colors[c(6,7,8)]
V(g_obs_S)$color <- adjustcolor(Customized_colors[RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_2d_R1_hat_z], alpha.f = .7)
par(mfrow=c(1,1),mai = c(0.05, 0.05, 0.05, 0.05),mgp = c(1,0.25,0))
plot(g_obs_S, rescale=T,layout=RDA_SampsonMonks_Directed_ZIPLPCM_Sup_Beta_1_9_T60k_2d_R1_hat_U,edge.curved=0.0,edge.width=E(g_obs_S)$weight*0.5,edge.arrow.size=0.1,vertex.frame.width=0.001)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# 2-dimensional visualization of Windsurfers
library("IMIFA")
U_ref <- RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_U
U_ref[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_z==1,] <- c(3,-3)
U_ref[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_z==2,] <- c(3,-3)
U_ref[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_z==3,] <- c(-3,3)
RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_U <- # rotate the whole plot to better match the pattern of 3-d visualization
  Procrustes(RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_U,
             U_ref, translate = TRUE ,dilate = FALSE)$X.new

My_colors <- c(brewer.pal(10,"RdBu")[c(4,8)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])
g_obs_W <- graph_from_adjacency_matrix(Windsurfers_adj,mode = "undirected",weighted = TRUE)
E(g_obs_W)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(Windsurfers_adj))[E(g_obs_W)$weight]
betw <- betweenness(g_obs_W)
V(g_obs_W)$size <- sqrt(betw/1.5+mean(betw))*0.8 # set the vertex size
V(g_obs_W)$frame.color <- "black"
V(g_obs_W)$label <- "" 
Customized_colors <- My_colors[c(1,6,2)]
V(g_obs_W)$color <- adjustcolor(Customized_colors[RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_z], alpha.f = .7)
par(mfrow=c(1,1),mai = c(0.05, 0.05, 0.05, 0.05),mgp = c(1,0.25,0))
plot(g_obs_W, rescale=T,layout=RDA_Windsurfers_UnDirected_ZIPLPCM_unSup_Beta_1_19_T60k_2d_R1_hat_U,edge.curved=0.0,edge.width=E(g_obs)$weight*0.25)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# 2-dimensional visualization of Train bombing suspects
g_obs_T <- graph_from_adjacency_matrix(Train_bombing_adj,mode = "undirected",weighted = TRUE)
E(g_obs_T)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(Train_bombing_adj))[E(g_obs_T)$weight]
betw <- betweenness(g_obs_T)
V(g_obs_T)$size <- sqrt(betw/1.5+mean(betw))*0.65
V(g_obs_T)$frame.color <- "black"
V(g_obs_T)$label <- "" 
Customized_colors <- My_colors[c(17)]
V(g_obs_T)$color <- adjustcolor(Customized_colors[RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_Beta_1_9_T60k_2d_R1_hat_z], alpha.f = .7)
par(mfrow=c(1,1),mai = c(0.05, 0.05, 0.05, 0.05),mgp = c(1,0.25,0))
plot(g_obs_T, rescale=T,layout=RDA_TrainBombing_UnDirected_ZIPLPCM_unSup_Beta_1_9_T60k_2d_R1_hat_U,edge.curved=0.0,edge.width=E(g_obs_T)$weight)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# 2-dimensional visualization of Summit co-attendance criminality network
library("igraph")
library("RColorBrewer")
My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[4],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])
g_obs_I <- graph_from_adjacency_matrix(RDA_criminalNet$Y,mode = "undirected",weighted = TRUE)
E(g_obs_I)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(RDA_criminalNet$Y))[E(g_obs_I)$weight]
betw <- betweenness(g_obs_I)
V(g_obs_I)$size <- sqrt(betw/1.5+mean(betw))*0.45
V(g_obs_I)$shape <- c("circle","square")[c(as.factor(RDA_criminalNet$role))]
V(g_obs_I)$frame.color <- "black"
V(g_obs_I)$label <- ""
par(mfrow=c(1,2),mai = c(0.05, 0.05, 0.05, 0.05),mgp = c(1,0.25,0))
# Make the plot of hat_U and z*
V(g_obs_I)$color <- adjustcolor(My_colors[c(as.factor(RDA_criminalNet$RoleAffiliation))], alpha.f = 0.7)
plot(g_obs_I, rescale=T,layout=-RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_Beta_1_6_T60k_DA_2d_R1_hat_U,edge.curved=0.0,edge.width=E(g_obs_I)$weight*0.25,vertex.frame.width=0.001)
# Make the plot of hat_U and hat_z
Customized_colors <- My_colors[c(3,6,9,4,7,5,8)]
V(g_obs_I)$color <- adjustcolor(Customized_colors[RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_Beta_1_6_T60k_DA_2d_R1_hat_z], alpha.f = .7)
plot(g_obs_I, rescale=T,layout=-RDA_CriminalSummit_UnDirected_ZIPLPCM_Sup_Beta_1_6_T60k_DA_2d_R1_hat_U,edge.curved=0.0,edge.width=E(g_obs_I)$weight*0.25,vertex.frame.width=0.001)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))
```

The **Figure 15** of the paper can thus be recovered.

Here we complete this tutorial of the coding for the **ZIP-LPCM** paper.



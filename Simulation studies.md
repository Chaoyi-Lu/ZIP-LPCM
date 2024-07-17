# Simulation Studies for ZIP-LPCM-MFM

This tutorial includes the coding for the two simulation studies illusrated in the paper "A Zero-Inflated Latent Position Cluster Model with Mixture of Finite Mixtures" (ZIP-LPCM-MFM).

Recall here that we have two different scenarios within each simulation study.
The scenario 1 in the first simulation study (SS1) focuses on a network randomly generated from a zero-inflated Poisson latent position cluster model (ZIP-LPCM) while the scenario 2 works on a network randomly generated from a Poisson latent position cluster model (Pois-LPCM).
In our second simulation study (SS2), we mainly focus on the networks randomly generated from the zero-inflated Poisson stochastic block model (ZIP-SBM), which is the one we newly proposed in [Lu, C., Durante, D., and Friel, N. [2024+], "Zero-inflated stochastic block modeling of efficiency-security tradeoffs in weighted criminal networks"]().
The synthetic network in scenario 2 of SS2 is equipped with a hub while the one in scenario 1 is not.

The source code of all the functions required for the experiments in our paper is included in the [`Functions.R`] file of this repository.
We load the functions via the code below.

``` r
rm(list=ls())
gc()
source("Functions.R")
```

Within the [`Functions.R`] file, we also added some comments to explain the corresponding lines of the code.
Here we refer to that file for more details.

## Simulation study 1

The scenario 1 network of this first simulation study is randomly generated from a ZIP-LPCM with the following code.

``` r
SS1_Scenario1_Directed_ZIPLPCM <-
  Simulation_Directed_ZIPLPCM(beta=3,P=matrix(c(0.4,0.1,0.05,0.1,0.05,  0.05,0.4,0.1,0.05,0.1,  0.1,0.05,0.4,0.1,0.05,
                                                0.05,0.1,0.05,0.4,0.1,  0.1,0.05,0.1,0.05,0.4),5,5),
                              mu=matrix(c(-1.5,-1.5,-1.5, -2,2,0, 2,-2,0, 2,2,-2, -2,-2,2),nrow=3,ncol=5),
                              tau=c(1/0.25,1/0.5,1/0.75,1/1,1/1.25),d=3,
                              z=c(rep(1,5),rep(2,10),rep(3,15),rep(4,20),rep(5,25)),seed=NULL)
```

where (i) "beta" corresponds to the intercept parameter $\beta$ of the model, (ii) "P" corresponds to the $\boldsymbol{P}$ which is a group-level $K \times K$ matrix indicating the probability of unusual zero for the interactions between each pair of groups, (iii) "mu" and "tau" correspond to the $\boldsymbol{\mu}$ and $\boldsymbol{\tau}$ which are, respectively, the group centres and the group precisions of the corresponding multivariate normal distributions, (iv) "z" corresponds to the latent clutering or membership variable $\boldsymbol{z}$ which is an $1\times N$ vector storing the pre-specified clustering. Here, $K$ is the number of non-empty clusters that can be easily extracted from $\boldsymbol{z}$.

The corresponding model parameters and latent variables are in line with those introduced in the ZIP-LPCM-MFM paper.
The above simulation function brings the following contents for a ZIP-LPCM network.

``` r
SS1_Scenario1_Directed_ZIPLPCM$Y
SS1_Scenario1_Directed_ZIPLPCM$nu
SS1_Scenario1_Directed_ZIPLPCM$z
SS1_Scenario1_Directed_ZIPLPCM$U
SS1_Scenario1_Directed_ZIPLPCM$Z
```

Here, "Y" corresponds to the $N \times N$ adjacency matrix $\boldsymbol{Y}$ of the network with each entry $y_{ij}$ indicating the interaction weight from node $i$ to node $j$; "nu" corresponds to the $N \times N$ unusual zero indicator matrix $\boldsymbol{\nu}$ with each entry $\nu_{ij}\in \\{0,1\\}$ indicating whether the corresponding $y_{ij}$ is an unusual zero ($\nu_{ij}=1$) or not ($\nu_{ij}=0$); "U" corresponds to the latent position variable $\boldsymbol{U}$ which is a $N \times d$ matrix with each row $i$ storing the corresponding latent position $\boldsymbol{u_i}$ assigned to the node $i$; "Z" corresponds to $\boldsymbol{Z}$ which is the $N\times K$ matrix form of the clustering $\boldsymbol{z}$.
Recall here that we treat the simulated latent variables and the above model parameters used for simulating the network as the references, denoted by $(\cdot)^*$.

These model parameters $\boldsymbol{\mu}^\*,\boldsymbol{\tau}^\*,\beta^\*$ and $\boldsymbol{P}^\*$ are also stored as variables in the code below.

``` r
SS1_Scenario1_Directed_ZIPLPCM_mu <- matrix(c(-1.5,-1.5,-1.5, -2,2,0, 2,-2,0, 2,2,-2, -2,-2,2),nrow=3,ncol=5)
SS1_Scenario1_Directed_ZIPLPCM_tau <- c(1/0.25,1/0.5,1/0.75,1/1,1/1.25)
SS1_Scenario1_Directed_ZIPLPCM_beta <- 3
SS1_Scenario1_Directed_ZIPLPCM_P <- matrix(c(0.4,0.1,0.05,0.1,0.05,  0.05,0.4,0.1,0.05,0.1,  0.1,0.05,0.4,0.1,0.05,  
                                             0.05,0.1,0.05,0.4,0.1,  0.1,0.05,0.1,0.05,0.4),5,5)
SS1_Scenario1_Directed_ZIPLPCM_p <- SS1_Scenario1_Directed_ZIPLPCM$Z%*%SS1_Scenario1_Directed_ZIPLPCM_P%*%t(SS1_Scenario1_Directed_ZIPLPCM$Z)
library(Rfast) # for Rfast::Dist()
SS1_Scenario1_Directed_ZIPLPCM_P_m0 <- (SS1_Scenario1_Directed_ZIPLPCM_p/(SS1_Scenario1_Directed_ZIPLPCM_p+(1-SS1_Scenario1_Directed_ZIPLPCM_p)*dpois(0,exp(SS1_Scenario1_Directed_ZIPLPCM_beta-Rfast::Dist(SS1_Scenario1_Directed_ZIPLPCM$U)))))*
  ((SS1_Scenario1_Directed_ZIPLPCM$Y+diag(1,nrow(SS1_Scenario1_Directed_ZIPLPCM$Y)))==0)
```

where `SS1_Scenario1_Directed_ZIPLPCM_p` above corresponds to an individual-level $N \times N$ matrix $\boldsymbol{p}$ with each entry being $p_{z_iz_j}$. 
While `SS1_Scenario1_Directed_ZIPLPCM_P_m0` is also an individual level $N \times N$ matrix with each entry corresponding to the statistic:

$$\text{P}(\nu_{ij}=1|y_{ij}=0,p_{z_iz_j},\beta,\boldsymbol{u_i},\boldsymbol{u_j})=\frac{p_{z_iz_j}}{p_{z_iz_j}+(1-p_{z_iz_j})f_{\text{Pois}}(0|\text{exp}(\beta-||\boldsymbol{u_i}-\boldsymbol{u_j}||)},$$ 

which is the conditional probability of an unusual zero provided that the corresponding observed interaction is a zero interaction.
Both of these two statistics are evaluated based on reference model parameters and latent variables here, and thus are treated as the corresponding reference statistics.

The simulated network as well as the simulated latent variables can be stored in separate `.csv` files following the code below.

``` r
write.csv(SS1_Scenario1_Directed_ZIPLPCM$Y,"Datasets/SS1_Scenario1_Directed_ZIPLPCM_Y.csv", row.names = FALSE)
write.csv(SS1_Scenario1_Directed_ZIPLPCM$nu,"Datasets/SS1_Scenario1_Directed_ZIPLPCM_nu.csv", row.names = FALSE)
write.csv(SS1_Scenario1_Directed_ZIPLPCM$z,"Datasets/SS1_Scenario1_Directed_ZIPLPCM_z.csv", row.names = FALSE)
write.csv(SS1_Scenario1_Directed_ZIPLPCM$U,"Datasets/SS1_Scenario1_Directed_ZIPLPCM_U.csv", row.names = FALSE)
write.csv(SS1_Scenario1_Directed_ZIPLPCM_A,"Datasets/SS1_Scenario1_Directed_ZIPLPCM_A.csv", row.names = FALSE)
```

where "A" here corresponds to the exogenous node attributes $\boldsymbol{c}$ we denote in the paper, and is a contaminated
version of the true clustering $\boldsymbol{z}^\*$, reallocating the clustering of 20 out of 75 nodes following the code:

``` r
SS1_Scenario1_Directed_ZIPLPCM_A <- SS1_Scenario1_Directed_ZIPLPCM$z
SS1_Scenario1_Directed_ZIPLPCM_A[sample(length(SS1_Scenario1_Directed_ZIPLPCM$z),20)] <- sample(1:5,20,replace=TRUE)
```

Once we stored the simulated network, we can reload the data again following the code:

``` r
SS1_Scenario1_Directed_ZIPLPCM <- 
  list(Y = as.matrix(read.csv("Datasets/SS1_Scenario1_Directed_ZIPLPCM_Y.csv",header = TRUE)),
       nu = as.matrix(read.csv("Datasets/SS1_Scenario1_Directed_ZIPLPCM_nu.csv",header = TRUE)),
       z = c(as.matrix(read.csv("Datasets/SS1_Scenario1_Directed_ZIPLPCM_z.csv",header = TRUE))),
       U = as.matrix(read.csv("Datasets/SS1_Scenario1_Directed_ZIPLPCM_U.csv",header = TRUE)),
       A = c(as.matrix(read.csv("Datasets/SS1_Scenario1_Directed_ZIPLPCM_A.csv",header = TRUE))))
SS1_Scenario1_Directed_ZIPLPCM$Z <- t(t(matrix(SS1_Scenario1_Directed_ZIPLPCM$z,length(SS1_Scenario1_Directed_ZIPLPCM$z),max(SS1_Scenario1_Directed_ZIPLPCM$z)))==(1:max(SS1_Scenario1_Directed_ZIPLPCM$z)))*1
colnames(SS1_Scenario1_Directed_ZIPLPCM$Y) <- c()
colnames(SS1_Scenario1_Directed_ZIPLPCM$nu) <- c()
colnames(SS1_Scenario1_Directed_ZIPLPCM$U) <- c()
```

We can visualize the network via the 3-dimensional interactive plot of the reference latent positions $\boldsymbol{U}^\*$, the reference clustering $\boldsymbol{z}^\*$ and those interactions stored in $\boldsymbol{Y}$.
The plot is uploaded on Github at [`/Interactive 3-d latent positions plots/SS1_Scenario1_InteractivePlot.html`], and can be reproduced following the code below.
However, the Github cannot directly load the [`SS1_Scenario1_InteractivePlot.html`] file on the page, so the readers need to download it first in order to play with the interactive 3-d plot.

``` r
library("igraph")
library("RColorBrewer")

My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

g_obs <- graph_from_adjacency_matrix(SS1_Scenario1_Directed_ZIPLPCM$Y,mode = "directed",weighted = TRUE)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(SS1_Scenario1_Directed_ZIPLPCM$Y))[E(g_obs)$weight]
betw <- betweenness(g_obs) # evaluate the betweeness
VertexSize <- sqrt(betw/1.5+mean(betw))*1 # set the vertex size

library("plotly")
fig <- plot_ly() %>% # plot the reference clustering and U
  add_markers(x = SS1_Scenario1_Directed_ZIPLPCM$U[,1],
              y = SS1_Scenario1_Directed_ZIPLPCM$U[,2],
              z = SS1_Scenario1_Directed_ZIPLPCM$U[,3],
              text=paste("Node:",1:nrow(SS1_Scenario1_Directed_ZIPLPCM$Y),"<br>z*:",SS1_Scenario1_Directed_ZIPLPCM$z),
              size=VertexSize,sizes=c(100,300),
              color=as.factor(SS1_Scenario1_Directed_ZIPLPCM$z),colors=My_colors[6:10]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],1],
              y = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],2],
              z = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.35*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "ZIP-LPCM U* and z*",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

Running the code above would bring some warning messages and these are raised by some package bugs which are not yet resolved.
However, these bugs would not affect our output and the readers can ignore those warning messages.
The above code for the interactive 3-d plot also attached some texts or notes for each node and each non-zero interaction to help readers have a better understanding of the network.
If the readers put the mouse pointer on each node of the interactive plot, there will be a comment bracket showing (i) the coordinate of the node, (ii) the node number (e.g. node 1, node 2, ...), (iii) the reference clustering of the node.
If the mouse pointer is put on each non-zero interaction, the bracket will show (i) either the start coordinate or the end coordinate of the interaction vector, (ii) the interaction weight, (iii) an indicator of whether the interaction is from node $i$ to node $j$ where $i\< j$, i.e., whether the corresponding $y_{ij}$ is the upper-diagonal entry of the $\boldsymbol{Y}$.

The Figure 1 in the ZIP-LPCM-MFM paper can be recovered following the code below where the functions `subplot()` and `orca()` are leveraged to aggregate multiple interactive plots from different specific angles in one figure and to output a high-quality screenshot, respectively.

``` r
g_obs <- graph_from_adjacency_matrix(SS1_Scenario1_Directed_ZIPLPCM$Y,mode = "directed",weighted = TRUE)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(SS1_Scenario1_Directed_ZIPLPCM$Y))[E(g_obs)$weight]
betw <- betweenness(g_obs)
VertexSize <- sqrt(betw/1.5+mean(betw))*1

library("plotly")
# First interactive plot
fig1 <- plot_ly(scene ="scene1") %>% 
  add_markers(x = SS1_Scenario1_Directed_ZIPLPCM$U[,1],
              y = SS1_Scenario1_Directed_ZIPLPCM$U[,2],
              z = SS1_Scenario1_Directed_ZIPLPCM$U[,3],
              text=paste("Node:",1:nrow(SS1_Scenario1_Directed_ZIPLPCM$Y),"<br>z*:",SS1_Scenario1_Directed_ZIPLPCM$z),
              size=VertexSize,sizes=c(100,300),showlegend = TRUE,
              color=as.factor(SS1_Scenario1_Directed_ZIPLPCM$z),colors=My_colors[6:10]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig1 <- fig1 %>%
    add_trace(x = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],1],
              y = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],2],
              z = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.35*E(g_obs)$weight[i]))
}
# Second interactive plot
fig2 <- plot_ly(scene ="scene2") %>% 
  add_markers(x = SS1_Scenario1_Directed_ZIPLPCM$U[,1],
              y = SS1_Scenario1_Directed_ZIPLPCM$U[,2],
              z = SS1_Scenario1_Directed_ZIPLPCM$U[,3],
              text=paste("Node:",1:nrow(SS1_Scenario1_Directed_ZIPLPCM$Y),"<br>z*:",SS1_Scenario1_Directed_ZIPLPCM$z),
              size=VertexSize,sizes=c(100,300),showlegend = FALSE,
              color=as.factor(SS1_Scenario1_Directed_ZIPLPCM$z),colors=My_colors[6:10]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig2 <- fig2 %>%
    add_trace(x = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],1],
              y = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],2],
              z = SS1_Scenario1_Directed_ZIPLPCM$U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.35*E(g_obs)$weight[i]))
}
# Combine two plots in one and set the default visualizing angles
fig <- subplot(fig1, fig2) 
fig <- fig %>% layout(title = "", margin = list(l = 0,r = 0,b = 0,t = 0,pad = 0),
                      scene = list(domain=list(x=c(0,1/2),y=c(0,1)),
                                   xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                   camera = list(eye = list(x = 1.25, y = 1.25, z = 1.25)), # set the default visualizing angle
                                   aspectmode='auto'),
                      scene2 = list(domain=list(x=c(1/2,0.999),y=c(0,1)),
                                    xaxis = list(title = ''),yaxis = list(title = ''),zaxis = list(title = ''),
                                    camera = list(eye = list(x = -1.25, y = 1.25, z = 1.25)),
                                    aspectmode='auto'))
fig
```

Running the `fig` above will directly output an interactive plot which contains two interactive subplots.
The two interactive subplots are exactly the same but with different default visualizing angles.
Direct screenshot of the interactive plot would bring an low-quality figure.
Thus we leverage the function `orca()` shown below to obtain a high-quality screenshot of the interactive plots.
However, before running the `orca()` function in the package `plotly`, the readers need to first install the "orca" app following the Github instructions: [https://plotly.com/r/static-image-export/](https://plotly.com/r/static-image-export/) and [https://github.com/plotly/orca#installation](https://github.com/plotly/orca#installation).

``` r
orca(fig, "SS1_Sce1_Obs.pdf",scale=1,width=1800,height=850)
```

The above code will output exactly the same figure as the Figure 1 of the ZIP-LPCM-MFM paper.

The following code provides some extra statistics about the interaction weights of the scenario 1 network.

``` r
# Number/proportion of true 0
sum(SS1_Scenario1_Directed_ZIPLPCM$nu[SS1_Scenario1_Directed_ZIPLPCM$Y == 0]==0)-75
(sum(SS1_Scenario1_Directed_ZIPLPCM$nu[SS1_Scenario1_Directed_ZIPLPCM$Y == 0]==0)-75)/(75*74)
# Number/proportion of unusual 0
sum(SS1_Scenario1_Directed_ZIPLPCM$nu[SS1_Scenario1_Directed_ZIPLPCM$Y == 0]==1)
(sum(SS1_Scenario1_Directed_ZIPLPCM$nu[SS1_Scenario1_Directed_ZIPLPCM$Y == 0]==1))/(75*74)
# Number/proportion of non 0
sum(SS1_Scenario1_Directed_ZIPLPCM$nu[SS1_Scenario1_Directed_ZIPLPCM$Y != 0]==0)
(sum(SS1_Scenario1_Directed_ZIPLPCM$nu[SS1_Scenario1_Directed_ZIPLPCM$Y != 0]==0))/(75*74)
# Distribution of interaction weights
table(SS1_Scenario1_Directed_ZIPLPCM$Y[!is.nan((SS1_Scenario1_Directed_ZIPLPCM$Y+diag(NaN,75)))])
# Plot the distribution of interaction weights
hist(SS1_Scenario1_Directed_ZIPLPCM$Y[!is.nan((SS1_Scenario1_Directed_ZIPLPCM$Y+diag(NaN,75)))],200,xlab = "",ylab = "", main = "")
```

We follow the similar simulation process above to simulate the scenario 2 network.
The scenario 2 network is randomly generated from a Pois-LPCM with the following code.

``` r
SS1_Scenario2_Directed_PoisLPCM <-
  Simulation_Directed_PoissonLPCM(beta=3,mu=matrix(c(-1.5,-1.5,-1.5, -2,2,0, 2,-2,0, 2,2,-2, -2,-2,2),nrow=3,ncol=5),
tau=c(1/0.25,1/0.5,1/0.75,1/1,1/1.25),d=3,z=c(rep(1,5),rep(2,10),rep(3,15),rep(4,20),rep(5,25)),seed=NULL)
```

This brings the following contents for a Pois-LPCM network.

``` r
SS1_Scenario1_Directed_PoisLPCM$Y
SS1_Scenario1_Directed_PoisLPCM$z
SS1_Scenario1_Directed_PoisLPCM$U
SS1_Scenario1_Directed_PoisLPCM$Z
```

These model parameters $\boldsymbol{\mu}^\*,\boldsymbol{\tau}^\*$ and $\beta^\*$ are stored as variables in the code below.

``` r
SS1_Scenario2_Directed_PoisLPCM_mu <- matrix(c(-1.5,-1.5,-1.5, -2,2,0, 2,-2,0, 2,2,-2, -2,-2,2),nrow=3,ncol=5)
SS1_Scenario2_Directed_PoisLPCM_tau <- c(1/0.25,1/0.5,1/0.75,1/1,1/1.25)
SS1_Scenario2_Directed_PoisLPCM_beta <- 3
```

The simulated network as well as the simulated latent variables can be stored in `.csv` files by the following code.

``` r
write.csv(SS1_Scenario2_Directed_PoisLPCM$Y,"Datasets/SS1_Scenario2_Directed_PoisLPCM_Y.csv", row.names = FALSE)
write.csv(SS1_Scenario2_Directed_PoisLPCM$z,"Datasets/SS1_Scenario2_Directed_PoisLPCM_z.csv", row.names = FALSE)
write.csv(SS1_Scenario2_Directed_PoisLPCM$U,"Datasets/SS1_Scenario2_Directed_PoisLPCM_U.csv", row.names = FALSE)
write.csv(SS1_Scenario2_Directed_PoisLPCM_A,"Datasets/SS1_Scenario2_Directed_PoisLPCM_A.csv", row.names = FALSE)
```

where we follow the similar manner as the scenario 1 to contaminate the true clustering to obtain the exogenous node attributes below.
The practitioners can also consider using the same node attributes as the scenario 1 here.

``` r
SS1_Scenario2_Directed_PoisLPCM_A <- SS1_Scenario2_Directed_PoisLPCM$z
SS1_Scenario2_Directed_PoisLPCM_A[sample(length(SS1_Scenario2_Directed_PoisLPCM$z),20)] <- sample(1:5,20,replace=TRUE)
```

Thus the data can be reloaded again via:

``` r
SS1_Scenario2_Directed_PoisLPCM <- 
  list(Y = as.matrix(read.csv("Datasets/SS1_Scenario2_Directed_PoisLPCM_Y.csv",header = TRUE)),
       z = c(as.matrix(read.csv("Datasets/SS1_Scenario2_Directed_PoisLPCM_z.csv",header = TRUE))),
       U = as.matrix(read.csv("Datasets/SS1_Scenario2_Directed_PoisLPCM_U.csv",header = TRUE)),
       A = c(as.matrix(read.csv("Datasets/SS1_Scenario2_Directed_PoisLPCM_A.csv",header = TRUE))))
SS1_Scenario2_Directed_PoisLPCM$Z <- t(t(matrix(SS1_Scenario2_Directed_PoisLPCM$z,length(SS1_Scenario2_Directed_PoisLPCM$z),max(SS1_Scenario2_Directed_PoisLPCM$z)))==(1:max(SS1_Scenario2_Directed_PoisLPCM$z)))*1
colnames(SS1_Scenario2_Directed_PoisLPCM$Y) <- c()
colnames(SS1_Scenario2_Directed_PoisLPCM$U) <- c()
```

The simulated network can be visualized via the 3-dimensional interactive plot uploaded on Github at [`/Interactive 3-d latent positions plots/SS1_Scenario2_InteractivePlot.html`], and can be reproduced following the code below.

``` r
library("igraph")
library("RColorBrewer")
My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])
g_obs <- graph_from_adjacency_matrix(SS1_Scenario2_Directed_PoisLPCM$Y,mode = "directed",weighted = TRUE)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(SS1_Scenario2_Directed_PoisLPCM$Y))[E(g_obs)$weight]
betw <- betweenness(g_obs) 
VertexSize <- sqrt(betw/1.5+mean(betw))*1
library("plotly")
fig <- plot_ly() %>% # plot the reference clustering and U
  add_markers(x = SS1_Scenario2_Directed_PoisLPCM$U[,1],
              y = SS1_Scenario2_Directed_PoisLPCM$U[,2],
              z = SS1_Scenario2_Directed_PoisLPCM$U[,3],
              text=paste("Node:",1:nrow(SS1_Scenario2_Directed_PoisLPCM$Y),"<br>z*:",SS1_Scenario2_Directed_PoisLPCM$z),
              size=VertexSize,sizes=c(100,300),
              color=as.factor(SS1_Scenario2_Directed_PoisLPCM$z),colors=My_colors[6:10]
  )
Edges <- get.edgelist(g_obs)
for (i in 1:nrow(Edges)){
  fig <- fig %>%
    add_trace(x = SS1_Scenario2_Directed_PoisLPCM$U[Edges[i,],1],
              y = SS1_Scenario2_Directed_PoisLPCM$U[Edges[i,],2],
              z = SS1_Scenario2_Directed_PoisLPCM$U[Edges[i,],3],
              text=paste("Weight:",E(g_obs)$weight[i],"<br>UpperDiag?:",Edges[i,1]<Edges[i,2]),
              type = "scatter3d", mode = "lines", showlegend = FALSE,line = list(color = E(g_obs)$color[i], width = 0.25*E(g_obs)$weight[i]))
}
fig <- fig %>% layout(title = "Pois-LPCM U* and z*",scene = list(xaxis = list(title = 'x1'),yaxis = list(title = 'x2'),zaxis = list(title = 'x3')))
fig
```

The adjacency matrices heatmap plots shown in Figure 2 of the ZIP-LPCM-MFM paper can be reproduced following the code:

``` r
library("RColorBrewer")
library("pheatmap")
My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"Reds")[c(9,6)],brewer.pal(9,"RdPu")[5],brewer.pal(9,"Greys")[c(3,6,9)],brewer.pal(9,"GnBu")[5])

SS1_Scenario1_Directed_ZIPLPCM_Y_dataframe <- as.data.frame(SS1_Scenario1_Directed_ZIPLPCM$Y)
rownames(SS1_Scenario1_Directed_ZIPLPCM_Y_dataframe) <- colnames(SS1_Scenario1_Directed_ZIPLPCM_Y_dataframe) <- 1:nrow(SS1_Scenario1_Directed_ZIPLPCM$Y)
annotation_row_z_ref <- as.data.frame(as.factor(matrix(SS1_Scenario1_Directed_ZIPLPCM$z,length(SS1_Scenario1_Directed_ZIPLPCM$z),1)))
colnames(annotation_row_z_ref) <- "z_ref"
annotation_colors_z_ref <- My_colors[c(6:10)]
names(annotation_colors_z_ref) <- sort(unique(annotation_row_z_ref$z_ref))
SS1_Scenario1_Directed_ZIPLPCM_Y_heatmap <-
  pheatmap(SS1_Scenario1_Directed_ZIPLPCM_Y_dataframe[order(SS1_Scenario1_Directed_ZIPLPCM$z),order(SS1_Scenario1_Directed_ZIPLPCM$z)],
           color=colorRampPalette(brewer.pal(9,"YlOrRd"))(max(SS1_Scenario1_Directed_ZIPLPCM$Y)+1),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(SS1_Scenario1_Directed_ZIPLPCM$z))!=0)),gaps_col=c(which(diff(sort(SS1_Scenario1_Directed_ZIPLPCM$z))!=0)))

SS1_Scenario2_Directed_PoisLPCM_Y_dataframe <- as.data.frame(SS1_Scenario2_Directed_PoisLPCM$Y)
rownames(SS1_Scenario2_Directed_PoisLPCM_Y_dataframe) <- colnames(SS1_Scenario2_Directed_PoisLPCM_Y_dataframe) <- 1:nrow(SS1_Scenario2_Directed_PoisLPCM$Y)
annotation_row_z_ref <- as.data.frame(as.factor(matrix(SS1_Scenario2_Directed_PoisLPCM$z,length(SS1_Scenario2_Directed_PoisLPCM$z),1)))
colnames(annotation_row_z_ref) <- "z_ref"
annotation_colors_z_ref <- My_colors[c(6:10)]
names(annotation_colors_z_ref) <- sort(unique(annotation_row_z_ref$z_ref))
SS1_Scenario2_Directed_PoisLPCM_Y_heatmap <-
  pheatmap(SS1_Scenario2_Directed_PoisLPCM_Y_dataframe[order(SS1_Scenario2_Directed_PoisLPCM$z),order(SS1_Scenario2_Directed_PoisLPCM$z)],
           color=colorRampPalette(brewer.pal(9,"YlOrRd"))(max(SS1_Scenario1_Directed_ZIPLPCM$Y)+1)[1:(max(SS1_Scenario2_Directed_PoisLPCM$Y)+1)],
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE,
           annotation_row = annotation_row_z_ref,annotation_col = annotation_row_z_ref,
           annotation_colors=list(z_ref=annotation_colors_z_ref),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_legend=FALSE,
           gaps_row=c(which(diff(sort(SS1_Scenario2_Directed_PoisLPCM$z))!=0)),gaps_col=c(which(diff(sort(SS1_Scenario2_Directed_PoisLPCM$z))!=0)))

library("grid")
library("gridExtra")
g <- grid.arrange(SS1_Scenario1_Directed_ZIPLPCM_Y_heatmap[[4]],
                  SS1_Scenario2_Directed_PoisLPCM_Y_heatmap[[4]],
                  nrow=1,ncol=2,vp=viewport(width=1, height=1))
Fig <- cowplot::ggdraw(g)
Fig
```

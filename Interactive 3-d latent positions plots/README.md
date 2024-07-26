In this [`Interactive 3-d latent positions plots`] file, 
we upload all the 3-dimensional interactive plots of the latent positions and clustering we illusrate 
in the the paper **"A Zero-Inflated Latent Position Cluster Model with Mixture of Finite Mixtures" (ZIP-LPCM-MFM)**.
Note that the GitHub page cannot directly have a view of these 3-d interactive plots which are stored as [`.html`] files, so the readers need to download these files first and then it should be straight-forward to see these 3-d plots by directly open these files by browser.

Within each 3-d interactive plot, the readers are free to rotate, zoom in or zoom out the plot.
There is also an operation panel placed on top-right of the 3-d interactive plot for readers to play with.
More details can be found in [https://plotly.com/r/3d-charts/](https://plotly.com/r/3d-charts/).

If the readers put the mouse pointer on each node of the 3-d interactive plot, there would be showing a comment bracket which contains some basic information about the corresponding node/individual including: (i) the coordinate of the node, (ii) the node number (e.g. node 1, node 2, ...).
Depending on the network data, some extra information might also be included, for example, the reference clustering, the exogenous node attributes and so on.

If the readers put the mouse pointer on each edge of the 3-d interactive plot, the comment bracket would show (i) either the start coordinate or the end coordinate of the interaction vector, (ii) the interaction weight.
If the network is directed, the comment bracket would also show a variable `UpperDiag?` which is a `TRUE` or `FALSE` variable indicating whether the corresponding $y_{ij}$ is placed at the upper-diagonal part of the adjacency matrix so that we could know the direction of such an edge.

Each 3-d interactive plot contained in this file, respectively, corresponds to:
<br>[`SS1_Scenario1_InteractivePlot.html`]: Simulation study 1 scenario 1 reference latent positions $`\boldsymbol{U}^*`$ and referernce clustering $`\boldsymbol{z}^*`$ shown as **Figure 1** of the **ZIP-LPCM-MFM** paper.
<br>[`SS1_Scenario2_InteractivePlot.html`]: Simulation study 1 scenario 2 reference latent positions $`\boldsymbol{U}^*`$ and referernce clustering $`\boldsymbol{z}^*`$ that are not shown in the paper.
<br>[`SS2_Scenario1_InteractivePlot.html`]: Simulation study 2 scenario 1 summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ shown as the 1st row plots of **Figure 5** in the paper.
<br>[`SS2_Scenario2_InteractivePlot.html`]: Simulation study 2 scenario 2 summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ shown as the 2nd row plots of **Figure 5** in the paper.


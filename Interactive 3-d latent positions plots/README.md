In this [`Interactive 3-d latent positions plots`] file, 
we upload all the 3-dimensional interactive plots of the latent positions and clustering we illusrate 
in the the paper **"A Zero-Inflated Poisson Latent Position Cluster Model" (ZIP-LPCM)**.
Note that the GitHub page cannot directly have a view of these 3-d interactive plots which are stored as [`.html`] files, so the readers need to download these files first and then it should be straight-forward to see these 3-d plots by directly open these files by browser.

Within each 3-d interactive plot, the readers are free to rotate, zoom in or zoom out the plot.
There is also an operation panel placed on top-right of the 3-d interactive plot for readers to play with.
More details can be found in [https://plotly.com/r/3d-charts/](https://plotly.com/r/3d-charts/).

If the readers put the mouse pointer on each node of the 3-d interactive plot, a comment bracket would appear within which some basic information about the corresponding node/individual is included: (i) the coordinate of the node, (ii) the node number (e.g. node 1, node 2, ...).
Depending on different network data, some extra information might also be included, for example, the reference clustering we ednote as $\boldsymbol{z}^*$, the exogenous node attributes $\boldsymbol{c}$ and so on.

If the readers put the mouse pointer on each edge of the 3-d interactive plot, the comment bracket would show (i) either the start coordinate or the end coordinate of the interaction vector, (ii) the interaction weight.
If the network is directed, the comment bracket would also show a variable `UpperDiag?` which is a `TRUE` or `FALSE` variable indicating whether or not the corresponding $y_{ij}$ is placed at the upper-diagonal part of the adjacency matrix $\boldsymbol{Y}$ for each edge, so that we could know the direction the edges.

Each 3-d interactive plot contained in this file, respectively, corresponds to:
<br>[`SS1_Scenario1_InteractivePlot.html`]: Simulation study 1 scenario 1 reference latent positions $`\boldsymbol{U}^*`$ and referernce clustering $`\boldsymbol{z}^*`$ shown as **Figure 1** of the **ZIP-LPCM** paper. 
<br>[`SS1_Scenario2_InteractivePlot.html`]: Simulation study 1 scenario 2 reference latent positions $`\boldsymbol{U}^*`$ and referernce clustering $`\boldsymbol{z}^*`$ that are not shown in the paper. 
<br>[`SS2_Scenario1_InteractivePlot.html`]: Simulation study 2 scenario 1 summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ shown as the 1st row plots of **Figure 5** in the paper. 
<br>[`SS2_Scenario2_InteractivePlot.html`]: Simulation study 2 scenario 2 summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ shown as the 2nd row plots of **Figure 5** in the paper.
<br>[`RDA_SampsonMonks_InteractivePlot.html`]: The summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ obtained by the real data application on the **Sampson Monks** directed real network shown as **Figure 7** of the paper.
<br>[`RDA_Windsurfers_InteractivePlot.html`]: The summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ obtained by the real data application on the **Windsurfers** undirected real network shown as **Figure 9** of the paper.
<br>[`RDA_TrainBombing_InteractivePlot.html`]: The summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ obtained by the real data application on the **Train Bombing** undirected real network shown as **Figure 11** of the paper.
<br>[`RDA_CriminalSummit_InteractivePlot_Ref_z.html`]: The summarized latent positions $`\hat{\boldsymbol{U}}`$ and reference clustering $`\boldsymbol{z}^*`$ in the real data application on the **'Ndrangheta Mafia** undirected real network shown as the 1st row plots in **Figure 13** of the paper.
<br>[`RDA_CriminalSummit_InteractivePlot_Hat_z.html`]: The summarized latent positions $`\hat{\boldsymbol{U}}`$ and summarized clustering $`\hat{\boldsymbol{z}}`$ in the real data application on the **'Ndrangheta Mafia** undirected real network shown as the 2nd row plots in **Figure 13** of the paper.


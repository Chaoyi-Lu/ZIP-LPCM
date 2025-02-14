# A Zero-Inflated Poisson Latent Position Cluster Model

This repository is a complementary material for the paper **"A Zero-Inflated Latent Position Cluster Model" (ZIP-LPCM)** to provide tutorials for the coding of the experiments illustrated in the paper.
This paper mainly works on the ZIP-LPCM incorporated with the **mixture-of-finite-mixtures (MFM)** clustering prior.
The inference is based on a novel **partially collapsed Metropolis-within-Gibbs (PCMwG)** sampler with a newly proposed **truncated absorb-eject (TAE)** move embedded inside.
A supervised version of the method is also provided.

The notebooks [`Simulation studies.md`] and [`Real data applications.md`] are the main focuses in this repository where step-by-step tutorials are provided.
In [`Simulation studies.md`], the code for the two simulation studies with two different scenarios within each of them illustrated in the **ZIP-LPCM** paper are provided.
While the [`Real data applications.md`] tutorial contains the code for getting access to those four publicly available real networks as well as the corresponding implementations and post-processing code.
All the required functions are included in the [`Functions.R`] file of this repository, and all the required data are in the [`Datasets`] file shown above.

The readers can directly get access to the 3-dimensional interactive plots for all the experiments in the [`Interactive 3-d latent positions plots`] file following the [`README.md`] file included inside.

# \src
This folder contains all scripts that were used in the analysis. 
Apart from the scripts that were used for each figure, the folder contains the following general scripts:

- [functions.R](functions.R) contains general functions that were used in the analysis

## Scripts that need to be run in order to reproduce the figures (in the following sequence)

- [getTFcheckpoint.R](getTFcheckpoint.R), [getMito.sh](getMito.sh) and [PrepareData.R](PrepareData.R) were used to convert the Cell Ranger output into an R object
- [QCAnalysis.R](QCAnalysis.R) runs the QCAnalysis and filters out low-quality cells
- [HVGandNorm.R](HVGandNorm.R) was used to identify HVGs and normalize the data
- [SNNCluster.R](SNNCluster.R), [SubCluster.R](SubCluster.R) and [CleanCluster.R](CleanCluster.R) (in this sequence) was used to perform the two-level clustering 
- [DEAnalysis.R](DEAnalysis.R) computes differential expression between all clusters for Supplementary Table 1
- [DiffusionMap.R](DiffusionMap.R) and [BranchDE.R](BranchDE.R) compute the diffusion map used for Figure 3 and 4 and test for differential expression along the branches
- [ProgenitorDE.R](ProgenitorDE.R) and [NPPIDE.R](NPPIDE.R) are used to compute differential expression between PI and NP progenitors and between NP/PI progenitors and the remaining luminal cells, respectively (Figure 5)
- [DiffusionMap_PI.R](DiffusionMap_PI.R) computes the embedding of PI cells on to the diffusion map (Figure 5)

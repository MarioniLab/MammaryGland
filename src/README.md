# \src
This folder contains all scripts that were used in the analysis. 
Apart from the scripts that were used for each figure, the folder contains the following general scripts:

- [functions.R](functions.R) contains general functions that were used in the analysis
- [getTFcheckpoint.R](getTFcheckpoint.R) and [PrepareData.R](PrepareData.R) were used to convert the Cell Ranger output into an R object
- [QCAnalysis.R](QCAnalysis.R), [Clustering.R](Clustering.R) and [ClusterBootstrap.R](ClusterBootstrap.R) (in this sequence) **need to be run before any of the figure scripts** as they produce necessary intermediate data
- [DEAnalysis.R](DEAnalysis.R) computes differential expression between all clusters for Supplementary Table 1

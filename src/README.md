# R code that was used in the analysis

- functions.R contains general functions that were used in the analysis
- getTFcheckpoint.R and PrepareData.R were used to convert the Cell Ranger output into an R object
- QCAnalysis.R, Clustering.R and ClusterBootstrap.R (in this sequence) need to be run before any of the figure scripts can be run as they produce necessary intermediate data
- DEAnalysis.R computes differential expression between all clusters for Supplementary Table 1

# Bootstrap for clustering

library(scran)
library(dplyr)
library(dynamicTreeCut)
require(doParallel)
require(fpc)
require(clValid)
source("functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# Gene and Cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]

# Normalization
m <- t(t(m)/pD$sf)

# Feature Selection
m.sub <- m[fD$highVar,]

# Combinations of clustering parameters
dms <- c("euclidean","pearson","spearman")
lks <- c("average","complete")
dss <- c(0,1,2)
nCores <- 3
cl <-makeCluster(nCores, type="FORK")
registerDoParallel(cl)
result <- foreach(i=seq_along(dss), .combine=c) %dopar% {
    ds <- dss[i]
    tmp0 <- list()
	for (lk in lks) {
	    for (dm in dms) {

		dis <- asDis(m.sub,dm=dm)
		res <- clusterboot(data=dis, B=100, clustermethod=dynamicCluster,
				   ds=ds, lk=lk)
		id <- paste(dm,lk,ds,sep="_")
		tmp0[[id]] <- res
	    }}
    return(tmp0)}
stopCluster(cl)
# Save
saveRDS(result, "../data/Robjects/ClusterBootstrap.rds")

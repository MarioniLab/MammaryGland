library(scran)
library(dplyr)
library(knitr)
library(ggplot2)
library(dynamicTreeCut)
library(Rtsne)
library(pheatmap)
source("functions.R")

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")

m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
#Gene and Cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]

#Normalization
clusters <- quickCluster(m)
pD$sf <- computeSumFactors(m,clusters=clusters)
m <- t(t(m)/pD$sf)
# 

#Run various clustering combinations

#Feature Selection
library(fpc)
m.sub <- t(m[fD$highVar,])
dms <- c("euclidean","pearson","spearman")
lks <- c("average","ward.D2","complete")
dss <- c(0,1,2)
require(doParallel)
require(fpc)
require(clValid)
nCores <- 3
cl <-makeCluster(nCores, type="FORK")
registerDoParallel(cl)
result <- foreach(i=seq_along(dss), .combine=c) %dopar% {
    ds <- dss[i]
    tmp0 <- list()
	for (lk in lks) {
	    for (dm in dms) {
		res <- clusterboot(data=m.sub, B=100, clustermethod=dynamicCluster,
				   dm=dm, ds=ds, lk=lk)
		id <- paste(dm,lk,ds,sep="_")
		tmp0[[id]] <- res
	    }}
    return(tmp0)}
stopCluster(cl)
saveRDS(result, "../data/Robjects/ClusterBootstrap.rds")

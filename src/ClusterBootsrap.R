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
fD <- mutate(fD, meanExpression=rowMeans(m),
	     vars=apply(m,1,var),
	     CV2=apply(m,1,cv2),
	     dm=DM(meanExpression,CV2))

brennecke <- BrenneckeHVG(m,fdr=0.1,minBiolDisp=0.25)

fD$brennecke <- fD$id %in% brennecke

#Run various clustering combinations

#Feature Selection
library(fpc)
m.sub <- m[fD$brennecke,]
dms <- c("euclidean","pearson")
lks <- c("average","ward.D2")
dss <- c(0,1)
require(doParallel)
require(fpc)
require(clValid)
nCores <- 2
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

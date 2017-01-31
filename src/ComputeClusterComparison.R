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
m <- t(t(m)/pD$sf)
# 

#Run various clustering combinations

fss <- c("highVar")
dms <- c("pearson","spearman","euclidean")
lks <- c("complete","average","ward.D2")
dss <- c(0,1,2,3,4)
require(doParallel)
require(fpc)
require(clValid)
nCores <- 3
cl <-makeCluster(nCores, type="FORK")
registerDoParallel(cl)
result <- foreach (i=seq_along(dms), .combine=rbind) %dopar% {
    dm <- dms[i]
    tmp0 <- NULL
    for (fs in fss) {
	for (lk in lks) {
	    for (ds in dss) {
		res <- compClustering(m,fD,fs=fs,dm=dm,lk=lk,ds=ds)
		tmp0 <- rbind(tmp0,res)
	    }}}
    return(tmp0)}
stopCluster(cl)

saveRDS(result,file="../data/Robjects/ClusterComparison.rds")

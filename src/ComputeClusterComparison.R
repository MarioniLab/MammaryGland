library(scran)
library(dplyr)
library(dynamicTreeCut)
source("functions.R")

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm.rds")

m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

#Gene and Cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]

#Normalization
m <- t(t(m)/pD$sf)
# 

#trafM
m <- t(log2(m[fD$highVar,]+1))

#Run various clustering combinations
gc()
fss <- c("highVar")
dms <- c("pearson","spearman","euclidean")
lks <- c("complete","average","ward.D2")
dss <- c(0,1,2,3,4)
require(doParallel)
require(fpc)
require(clValid)
nCores <- 6
cl <-makeCluster(nCores, type="FORK")
registerDoParallel(cl)
result <- foreach (i=seq_along(dms), .combine=rbind) %dopar% {
    dm <- dms[i]
    tmp0 <- NULL
    for (fs in fss) {
	for (lk in lks) {
	    for (ds in dss) {
		res <- compClustering(trafM=m,fs=fs,dm=dm,lk=lk,ds=ds)
		tmp0 <- rbind(tmp0,res)
	    }}}
    return(tmp0)}
stopCluster(cl)

saveRDS(result,file="../data/Robjects/secondRun_2500/ClusterComparison.rds")

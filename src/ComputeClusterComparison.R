library(scran)
library(dplyr)
library(dynamicTreeCut)
source("functions.R")

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm.rds")
set.seed(300)
dataList <- subSample(dataList,cell.number=10000)

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
dms <- c("pearson","spearman","euclidean")
lks <- c("average")
dss <- c(0,1,2,3,4)
require(doParallel)
require(fpc)
require(clValid)
nCores <- detectCores()
cl <-makeCluster(nCores, type="FORK")
registerDoParallel(cl)
result <- foreach (i=seq_along(dms)) %dopar% {
    dm <- dms[i]
    tmp0 <- list() 
	for (lk in lks) {
	    for (ds in dss) {
		res <- compClustering(trafM=m, dm=dm, lk=lk, ds=ds)
		tmp0[[paste(dm,lk,ds,sep="_")]] <- res
	    }}
    return(tmp0)}
stopCluster(cl)

result <- unlist(result,recursive=FALSE)
saveRDS(result,file="../data/Robjects/secondRun_2500/ClusterComparison.rds")

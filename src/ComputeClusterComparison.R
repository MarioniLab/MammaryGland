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

fit <- trendVar(log2(m+1))
decVar <- decomposeVar(log2(m+1),fit)
plot(decVar$total~decVar$mean)
feat <- decVar[decVar$FDR < 0.1,] %>% rownames()

highVar <- dplyr::arrange(fD,desc(dm)) %>% .$id
topX <- floor(0.05*length(highVar))

library(M3Drop)
brennecke <- BrenneckeGetVariableGenes(m,fdr=0.1,minBiolDisp=0.25)

fD$highDm <- fD$id %in% highVar[1:topX]
fD$highVar <- fD$id %in% feat
fD$brennecke <- fD$id %in% brennecke

#Run various clustering combinations

fss <- c("brennecke")
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

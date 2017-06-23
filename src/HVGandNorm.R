# Estimate size factors using scran
library(dplyr)
library(scran)

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
rm(dataList)

# Gene and cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]

clusters <- quickCluster(m,method="igraph")
pD$sf <- computeSumFactors(m,clusters=clusters)

# Normalize count matrix
m <- t(t(m)/pD$sf)


# Highly variable genes
brennecke <- BrenneckeHVG(m,fdr=0.1)
fD$highVar <- fD$id %in% brennecke

# Compute tSNE 
fPCA <- log2(t(m[brennecke,])+1)
fPCA <- scale(fPCA,scale=TRUE,center=TRUE)
set.seed(300)
tsn <- Rtsne(fPCA,perplexity=50)
pD$tSNE1 <- tsn$Y[,1]
pD$tSNE2 <- tsn$Y[,2]

# save
pD.add <- pD[,c("barcode","sf","tSNE1","tSNE2")]
fD.add <- fD[,c("id","highVar")]


dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC.rds")
pD <- left_join(dataList[[2]],pD.add)
fD <- left_join(dataList[[3]],fD.add)
saveRDS(dataList,file="../data/Robjects/secondRun_2500/test.rds")

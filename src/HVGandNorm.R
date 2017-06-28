# Estimate size factors using scran
library(dplyr)
library(scran)
library(Rtsne)
source("functions.R")

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC.rds")
# set.seed(300)
# dataList <- subSample(dataList, cell.number=2000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

# fD$keep <- rowMeans(m) > 0.01

# Gene and cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]

clusters <- quickCluster(m,method="hclust")
minSize <- min(table(clusters))
pD$sf <- computeSumFactors(m, sizes=seq(20,min(100,minSize),5),clusters=clusters)

plot(log10(colSums(m))~log10(pD$sf),main="Library Size versus Size Factors")

# Normalize count matrix
m <- t(t(m)/pD$sf)


# Highly variable genes 
varDf <- technicalCV2(t(t(m)/(1/pD$sf)), is.spike=NA, sf.cell=pD$sf, sf.spike=pD$sf)
plot(varDf$mean, varDf$cv2, log="xy",pch=19)
points(varDf$mean, varDf$trend, col="red", pch=16, cex=0.5)
points(varDf[varDf$FDR < 0.1,"mean"], varDf[varDf$FDR < 0.1, "cv2"], col="blue", pch=16, cex=0.5)

# Select only correlated genes
param <- MulticoreParam(workers=detectCores())
cors <- correlatePairs(m[rownames(varDf[varDf$FDR < 0.1,]),],BPPARAM=param,
		       per.gene=TRUE)
highVar.genes <- cors[cors$FDR <= 0.1, "gene"]

fD$highVar <- fD$id %in% highVar.genes



# Compute tSNE 
fPCA <- log2(t(m[fD$highVar,])+1)
fPCA <- scale(fPCA,scale=TRUE,center=TRUE)
set.seed(300)
tsn <- Rtsne(fPCA,perplexity=50)
pD$tSNE1 <- tsn$Y[,1]
pD$tSNE2 <- tsn$Y[,2]

# save
pD.add <- pD[,c("barcode","sf","tSNE1","tSNE2")]
fD.add <- fD[,c("id","highVar")]

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC.rds")
dataList[["phenoData"]] <- left_join(dataList[["phenoData"]],pD.add)
dataList[["featureData"]] <- left_join(dataList[["featureData"]],fD.add)

saveRDS(dataList,file="../data/Robjects/secondRun_2500/ExpressionList_QC_norm.rds")

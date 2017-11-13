# Estimate size factors using scran
library(dplyr)
library(scran)
library(Rtsne)
source("functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
# set.seed(300)
# dataList <- subSample(dataList, cell.number=5000)
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
var.des <- trendVar(log2(m+1),trend="semiloess")
var.out <- decomposeVar(log2(m+1),var.des)
o <- order(var.out$mean)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >=0.5),]
points(hvg.out$mean, hvg.out$total, pch=16, col="red")

# Select only correlated genes
# param <- MulticoreParam(workers=6)
# cors <- correlatePairs(m[rownames(var.out[var.out$FDR < 0.01,]),],BPPARAM=param,
#                        per.gene=TRUE)
# highVar.genes <- cors[cors$FDR <= 0.01, "gene"]

# add info to fD
stopifnot(identical(rownames(var.out),fD$id))
fD$highVar <- fD$id %in% rownames(hvg.out)
fD$highVarFDR <- var.out$FDR
fD$highVarBiolComp <- var.out$bio


# Compute tSNE 
fPCA <- log2(t(m[fD$highVar,])+1)
fPCA <- scale(fPCA,scale=TRUE,center=TRUE)
set.seed(300)
tsn <- Rtsne(fPCA,perplexity=50)
pD$tSNE1 <- tsn$Y[,1]
pD$tSNE2 <- tsn$Y[,2]

# save
pD.add <- pD[,c("barcode","sf","tSNE1","tSNE2")]
fD.add <- fD[,c("id","highVar","highVarFDR","highVarBiolComp")]

dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
dataList[["phenoData"]] <- left_join(dataList[["phenoData"]],pD.add)
fD.new <- left_join(dataList[["featureData"]],fD.add)
fD.new$highVar[is.na(fD.new$highVar)] <- FALSE
dataList[["featureData"]] <- fD.new

saveRDS(dataList,file="../data/Robjects/ExpressionList_QC_norm.rds")

# Estimate size factors using scran
library(dplyr)
library(scran)
library(Rtsne)
source("functions.R")

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered.rds")
# set.seed(300)
# dataList <- subSample(dataList, cell.number=2000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
# rm(dataList)

# Gene and cell filtering
m <- m[,pD$PassAll]
pD <- pD[pD$PassAll,]

# Normalize count matrix
m <- t(t(m)/pD$sf)

library(cluster)
# Define clustering function
func <- function(x,k) {
    dis <- as.dist((1-cor(t(x)))/2)
    tree <- hclust(dis, method="average")
    cluster <- cutree(tree,k=k)
    out <- list()
    out$cluster <- cluster
    return(out)
}

# sbst

# Perform subclustering
out <- NULL
for (cl in c(1:13)) {
    pD.sub <- pD[pD$Cluster==cl,]
    m.sub <- m[,pD.sub$barcode]
    m.sub <- m.sub[rowMeans(m.sub) >0.01,]

    # Highly variable genes
    var.des <- trendVar(log2(m.sub+1),trend="semiloess")
    var.out <- decomposeVar(log2(m.sub+1),var.des)
    o <- order(var.out$mean)
    hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >=0.5),]

    x <- log2(m.sub[rownames(hvg.out),]+1)
    gpas <- clusGap(t(x), func, K.max=3, B=100)
    k.opt <- maxSE(gpas$Tab[,"gap"],gpas$Tab[,"SE.sim"])
    subclust <- func(t(x), k.opt)
    pD.sub$SubCluster <- as.factor(paste(cl,subclust$cluster,sep="."))
    out <- rbind(out,pD.sub[,c("barcode","SubCluster")])
}

pD.add <- out

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered.rds")
dataList[["phenoData"]] <- left_join(dataList[["phenoData"]],pD.add)
saveRDS(dataList,file="../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered.rds")

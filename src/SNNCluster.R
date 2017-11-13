# Estimate size factors using scran
library(igraph)
library(scran)
library(dplyr)
source("functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

# Gene and cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]

# Normalize count matrix
m <- t(t(m)/pD$sf)

# SNN graph
igr <- buildSNNGraph(log2(m[fD$highVar,]+1),k=20)
cs <- cluster_louvain(igr)

pD$Cluster <- as.factor(cs$membership)

library(ggplot2)
library(cowplot)
p1 <- ggplot(pD, aes(tSNE1, tSNE2, color=Cluster)) +
    geom_point()

p2 <- ggplot(pD, aes(tSNE1, tSNE2, color=Cluster)) +
    geom_point() +
    facet_wrap(~Cluster)
plot_grid(p1,p2)

# save
pD.add <- pD[,c("barcode","Cluster")]
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm.rds")
dataList[["phenoData"]] <- left_join(dataList[["phenoData"]],pD.add)
saveRDS(dataList,file="../data/Robjects/ExpressionList_QC_norm_clustered.rds")

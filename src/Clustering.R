#Clustering
library(plyr)
library(scran)
library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
library(dynamicTreeCut)
library(Rtsne)
library(pheatmap)
source("functions.R")

#Reading Data
rnd_seed <- 300
clustDat <- readRDS("../data/Robjects/ClusterComparison.rds") 
bootDat <- readRDS("../data/Robjects/ClusterBootstrap.rds")
dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")

m.full <- dataList[[1]]
pD.full <- dataList[[2]]
fD.full <- dataList[[3]]

#Gene and Cell filtering
m <- m.full[fD.full$keep,pD.full$PassAll]
pD <- pD.full[pD.full$PassAll,]
fD <- fD.full[fD.full$keep,]

#Normalization
m <- t(t(m)/pD$sf)

#Clustering on highVar genes
m.sub <- m[fD$highVar,]

dm <- "euclidean"
lk <- "average"
ds <- 1
minSize <- 15

clust <- dynamicCluster(t(m.sub), dm=dm,
				 lk=lk,
				 ds=ds,
				 output="Cluster",
				 minSize=minSize)

pD$cluster <- clust$cluster

#Rename clusters according to nomenclature in paper
pD$cluster <- as.factor(mapvalues(pD$cluster,c(0,1,2,3,4,5,6,7,9,10,8),
			c(0,1,2,3,4,5,6,7,8,9,10)))

#Clusters per condition and sample
table(pD$cluster,pD$Condition)
table(pD$cluster,pD$SampleID)

# Removal of Immune Cells from the DataSet
pD.add <- select(pD, barcode, cluster)
pD.new <- left_join(pD.full,pD.add, by="barcode")
stopifnot(identical(colnames(m.full),
		    as.character(pD.full$barcode)))
stopifnot(identical(rownames(m.full),
		    as.character(fD.full$id)))
pD.new <- filter(pD.new, !(cluster %in% c(0,10)) & PassAll)
m.new <- m.full[,as.character(pD.new$barcode)]
fD.new <- fD.full

#Filter
isexpThreshold <- 10
expThreshold <- 50*isexpThreshold
keep1 <- rowSums(m.new!=0) > isexpThreshold
keep2 <- rowSums(m.new) > expThreshold
fD.new$keep <- keep1 & keep2
m.filtered <- m.new[fD.new$keep,]

#Normalize
clusters <- quickCluster(m.filtered)
pD.new$sf <- computeSumFactors(m.filtered,clusters=clusters)
plot(log10(colSums(m.filtered))~log10(pD.new$sf),main="Library Size versus Size Factors")
m.norm <- t(t(m.filtered)/pD.new$sf)

#HVG definition by Brennecke
brennecke <- BrenneckeHVG(m.norm,fdr=0.1)
fD.new$highVar <- fD.new$id %in% brennecke

#Dim Reduction
fPCA <- log2(t(m.norm[brennecke,])+1)
fPCA <- scale(fPCA,scale=TRUE,center=TRUE)

set.seed(rnd_seed)
tsn <- Rtsne(fPCA,perplexity=25)

pD.new$tSNE1 <- tsn$Y[,1]
pD.new$tSNE2 <- tsn$Y[,2]

# Save Data
pD.full <- left_join(pD.full, select(pD, cluster, barcode))
pD.addTo <- select(pD.full, -sf) %>% 
    mutate(isImmuneCell=cluster==10,
	   isOutlier=cluster==0)
pD.add <- select(pD.new, barcode, sf, tSNE1, tSNE2)
pD.out <- left_join(pD.addTo,pD.add, by="barcode")

fD.addTo <- select(fD.full, -highVar, -keep)
fD.add <- select(fD.new, highVar, keep, id)
fD.out <- left_join(fD.addTo,fD.add, by="id")

stopifnot(identical(colnames(m.full),
		    as.character(pD.out$barcode)))
stopifnot(identical(rownames(m.full),
		    as.character(fD.out$id)))
out <- list()
out[[1]] <- m.full
out[[2]] <- pD.out
out[[3]] <- fD.out
saveRDS(out,file="../data/Robjects/ExpressionList_Clustered.rds")

# DE for all clusters to find marker genes

library(scran)

# Load Data
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]

# Cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

# norm
m <- t(t(m)/pD$sf)

# Genes
grps <- unique(pD$SuperCluster)
keep <- NULL
for (grp in grps) {
    sbst <- as.character(pD[pD$SuperCluster==grp,"barcode"])
    n <- m[,sbst]
    tmp <- rownames(n[rowMedians(n)>=1,])
    keep <- c(tmp,keep)
}
keep <- rownames(m) %in% unique(keep)
m <- m[keep,]
fD <- fD[keep,]

# log
m <- log2(m+1)

# markers
cls <- as.character(pD$SuperCluster)
markers <- findMarkers(m,cls)

library(dplyr)
geneIDs <- fD[,c("id","symbol")]
colnames(geneIDs) <- c("Gene","Symbol")
out <- lapply(markers, function(d) {
		    out <- left_join(d,geneIDs)
		    rownames(out) <- out$Gene
		    return(out)
})


saveRDS(markers,"../data/Robjects/secondRun_2500/DEList.rds")

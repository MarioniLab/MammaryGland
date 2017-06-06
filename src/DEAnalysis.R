# DE for all clusters to find marker genes

library(scran)

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# Cells
keepCells <- pD$PassAll & !(pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

# norm
m <- t(t(m)/pD$sf)

# Genes
grps <- unique(pD$cluster)
keep <- NULL
for (grp in grps) {
    sbst <- as.character(pD[pD$cluster==grp,"barcode"])
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
cls <- as.character(pD$cluster)
markers <- findMarkers(m,cls)

saveRDS(markers,"../data/Robjects/DEList.rds")

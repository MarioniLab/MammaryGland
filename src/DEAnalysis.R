# DE for all clusters to find marker genes

library(scran)

# Load Data
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]

# Cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

# norm
m <- t(t(m)/pD$sf)

# Genes
grps <- unique(pD$SubClusterNumbers)
keep <- NULL
for (grp in grps) {
    sbst <- as.character(pD[pD$SubClusterNumbers==grp,"barcode"])
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
cls <- as.character(pD$SubClusterNumbers)
markers <- findMarkers(m,cls)

library(dplyr)
geneIDs <- fD[,c("id","symbol")]
colnames(geneIDs) <- c("Gene","Symbol")
out <- lapply(markers, function(d) {
		    out <- left_join(d,geneIDs)
		    rownames(out) <- out$Gene
		    colnames(out) <- gsub("^C","logFC.vs.C",colnames(out))
		    return(out)
})

options(java.parameters="-Xmx16000m")
library(xlsx)
jgc <- function()
{
  .jcall("java/lang/System", method = "gc")
} 

for (cluster in names(out)) {
    write.xlsx(data.frame(out[[cluster]]), file="../data/MarkerGenes.xlsx", sheetName=cluster, row.names=FALSE, append = TRUE)
    jgc()
}

saveRDS(markers,"../data/Robjects/DEList.rds")

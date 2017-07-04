# DE for all clusters to find marker genes

library(scran)
library(Rtsne)
library(limma)
source("functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered2.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
rm(dataList)

# Cells
keepCells <- pD$PassAll 
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
m <- m[fD$keep,]
fD <- fD[fD$keep,]


# Norm
m <- t(t(m)/pD$sf)

# Log
m <- log2(m+1)

# Markers
cls <- as.character(pD$SubCluster)
markers <- deAll(m,cls)

out <- NULL
for (mn in names(markers)) {
    deList <- markers[[mn]]
    tmp0 <- apply(deList[,grepl("FDR",colnames(deList))],2, function(x) sum(x<0.01))
    tmp1 <- gsub("FDR.","",names(which(tmp0 <= 10)))
    if (length(tmp1>0)){ names(tmp1) <- mn}
    out <- c(out,tmp1)
}

print(out)


#Renaming clusters and merging 7.1 and 8.1

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered2.rds")
pD <- dataList[[2]]

library(plyr)
pD$SubCluster <- as.character(pD$SubCluster)
pD$SubCluster <- mapvalues(pD$SubCluster, 
			   c("1.1","2.1","3.1","3.2","4.1","5.1","5.2",
			     "6.1","6.2","7.1","8.1","9.1","9.2","10.1",
			     "10.2","10.3","11.1","11.2","12.1","12.2","13.1"),
			   c("C5","C1-PI","C3-PI","C3-NP","C4","C1-NP","C1-L",
			     "6-1","6-2","C7","C7","C8early","C8late","10.1","C9",
			     "10.3","C6-G1","C6-G2","C2diff","C2pro","C6"))

pD$SubCluster <- factor(pD$SubCluster, levels=c("C1-NP","C1-PI","C1-L","C2diff","C2pro",
						"C3-NP","C3-PI","C4","C5","C6",
						"C6-G1","C6-G2","C7","C8early","C8late",
						"C9","6-1","6-2","10.1","10.3"))

dataList[[2]] <- pD

#############Two scripts merged here needs polishing

#reload data
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
rm(dataList)

# Mark cells
pD$IsNonEpithelial <- pD$SubCluster %in% c("6-1","6-2","10.1","10.3")

# Gene and cell filtering
m <- m[,pD$PassAll & !pD$IsNonEpithelial]
pD <- pD[pD$PassAll & !pD$IsNonEpithelial,]
fD$keep <- rowMeans(m)>0.01
m <- m[fD$keep,]
fD <- fD[fD$keep,]

# ---- RecomputeHVGs ----

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
library(dplyr)
xchange.pD <- c("tSNE1","tSNE2")
xchange.fD <- c("keep","highVar","highVarFDR","highVarBiolComp")
pD.add <- pD[,c("barcode",xchange.pD)]
fD.add <- fD[,c("id", xchange.fD)]

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered2.rds")
pD.old <- dataList[["phenoData"]]
pD.old <- pD.old[,!(colnames(pD.old) %in% xchange.pD)]
pD.new <- left_join(pD.old,pD.add)
pD.new$IsNonEpithelial <- pD.new$SubCluster %in% c("6-1","6-2","10.1","10.3")
pD.new$IsNonEpithelial[is.na(pD.new$IsNonEpithelial)] <- FALSE
dataList[["phenoData"]] <- pD.new

fD.old <- dataList[["featureData"]]
fD.old <- fD.old[,!(colnames(fD.old) %in% xchange.fD)]
fD.new <- left_join(fD.old,fD.add)
fD.new$highVar[is.na(fD.new$highVar)] <- FALSE
fD.new$keep[is.na(fD.new$keep)] <- FALSE
dataList[["featureData"]] <- fD.new

saveRDS(dataList,file="../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered3.rds")

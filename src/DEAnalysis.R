library(scran)
library(edgeR)
library(dplyr)
library(ggplot2)
library(Rtsne)
library(doParallel)
source("functions.R")

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

# Genes
keep <- rowMeans(m) > 0.1
m <- m[keep,]
fD <- fD[keep,]


# DGEList
nf <- log(pD$sf/pD$UmiSums)
pD$nf <- exp(nf-mean(nf))
y <- DGEList(counts=m,
	     samples=pD,
	     genes=fD,
	     norm.factors=pD$nf)

# Model and disp
cluster <- factor(pD$cluster)
de.design <- model.matrix(~0+cluster)
y <- estimateDisp(y, de.design, trend="none", df=0)
fit <- glmFit(y, de.design)

nCores <- 3
cl <-makeCluster(nCores, type="FORK")
registerDoParallel(cl)
result <- foreach (i=seq_along(levels(cluster))) %dopar% {
    result.logFC <- result.PValue <- list()
    chosen.clust <- i

    for ( clust in seq_len(nlevels(cluster))) {
	if (clust==chosen.clust) { next } 
	print(paste0("We are at cluster ",clust))
	contrast <- numeric(ncol(de.design))
	contrast[chosen.clust] <- 1
	contrast[clust] <- -1
	res <- glmTreat(fit, contrast=contrast, lfc=1)
	con.name <- paste0('vs.', levels(cluster)[clust])
	result.logFC[[con.name]] <- res$table$logFC
	result.PValue[[con.name]] <- rank(res$table$PValue, ties="first")
    }

    collected.ranks <- lapply(result.PValue, rank, ties="first")
    min.rank <- do.call(pmin, collected.ranks)
    marker.set <- data.frame(Top=min.rank, Gene=rownames(y),
			     logFC=do.call(cbind, result.logFC),
			     top=do.call(cbind, result.PValue),
			     stringsAsFactors=FALSE)
    marker.set <- marker.set[order(marker.set$Top),]
    return(marker.set)
}
names(result) <- levels(cluster)

saveRDS(result,"../data/Robjects/DEList.rds")

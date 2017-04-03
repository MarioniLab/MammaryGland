library(scran)
library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
library(dynamicTreeCut)
library(Rtsne)
library(pheatmap)
source("functions.R")

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- pD$PassAll & !(pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep <- rowMeans(m) > 0.1
m <- m[keep,]
fD <- fD[keep,]


library(edgeR)
nf <- log(pD$sf/pD$UmiSums)
pD$nf <- exp(nf-mean(nf))
y <- DGEList(counts=m,
	     samples=pD,
	     genes=fD,
	     norm.factors=pD$nf)

cluster <- factor(pD$cluster)
de.design <- model.matrix(~0+cluster)
y <- estimateDisp(y, de.design)
fit <- glmFit(y, de.design)

#copy cat
library(doParallel)
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
names(result) <- seq_along(levels(cluster))

saveRDS(result,"../data/Robjects/DEList.rds")
#Select only up regulated
# marker.set <- result[[2]]
# sbst <- marker.set[,-c(1,2)]
# allUp <- apply(sbst,1, function(x) sum(x > 0)==length(x))
# marker.up <- marker.set[allUp,]
# marker.up <- marker.up[order(apply(marker.up[,-c(1,2)],1,min),decreasing=TRUE),]
# 
# m.n <- t(t(m)/pD$sf)
# gns <- x$Gene[1:128]
# plotGeneDist(m.n,pD,fD,gns,colorBy="cluster")
# x <- marker.set[order(abs(marker.set$logFC.vs.4),decreasing=TRUE),]
# x
# 
#SCDE test
# library(scde)
# m <- apply(m,2, function(x) {storage.mode(x) <- 'integer';x})
# o.ifm <- scde::scde.error.models(
#     counts = m,
#     groups = cluster,
#     n.cores = 3,
#     threshold.segmentation = TRUE,
#     save.crossfit.plots = FALSE,
#     min.count.threshold= 1,
#     save.model.plots = FALSE,
#     verbose = 1,
#     min.size.entries = 2
# )
# 
# priors <- scde::scde.expression.prior(
#     models = o.ifm,
#     counts = m,
#     length.out = 400,
#     show.plot = FALSE
# )
# 
# resSCDE <- scde::scde.expression.difference(
#     o.ifm,
#     cnts,
#     priors,
#     groups = cond,
#     n.randomizations = 100,
#     n.cores = 1,
#     verbose = 0
# )

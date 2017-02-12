library(scran)
library(scde)
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
keepCells <- pD$PassAll & pD$cluster !=0
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep1 <- rowSums(m!=0)>10
keep2 <- rowSums(m)>20
keep <- keep1 & keep2
m <- m[keep,]
fD <- fD[keep,]

#SCDE test
library(scde)
m <- apply(m,2, function(x) {storage.mode(x) <- 'integer';x})
cluster <- factor(pD$cluster)
o.ifm <- scde::scde.error.models(
				 counts = m,
				 groups = cluster,
				 n.cores = 1,
				 threshold.segmentation = TRUE,
				 save.crossfit.plots = FALSE,
				 min.count.threshold= 1,
				 save.model.plots = FALSE,
				 verbose = 1,
				 min.size.entries = 2
				 )

priors <- scde::scde.expression.prior(
				      models = o.ifm,
				      counts = m,
				      length.out = 400,
				      show.plot = FALSE
				      )

resSCDE <- scde::scde.expression.difference(
					    o.ifm,
					    cnts,
					    priors,
					    groups = cond,
					    n.randomizations = 100,
					    n.cores = 1,
					    verbose = 0
					    )

saveRDS(result,"../data/Robjects/DEList.rds")

# Differentation inference with monocle

library(plyr)
library(dplyr)
library(ggplot2)
library(monocle)
library(RColorBrewer)
library(scran)
source("functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# ---- Prepocessing ----

# Cells
keepCells <- pD$keep & pD$Condition %in% c("NP","G") & !(pD$SuperCluster %in% c("C6","C6-G1","C7","C9"))
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep <- rowMeans(m)>0.01
m <- m[keep,]
fD <- fD[keep,]

# Normalize
m.norm <- t(t(m)/pD$sf)
rownames(m.norm) <- fD$symbol

# High Var
var.des <- trendVar(log2(m.norm+1),trend="semiloess")
var.out <- decomposeVar(log2(m.norm+1),var.des)
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >=0.5),]
fD$highVar <- fD$id %in% rownames(hvg.out)

# ---- Monocle ----

# Initiate Monocle
pD.af <- new("AnnotatedDataFrame", data=pD)
sampleNames(pD.af) <- pD$barcode
fD.af <- new("AnnotatedDataFrame", data=fD)
sampleNames(fD.af) <- fD$id

cds <- newCellDataSet(as.matrix(m),
	       phenoData=pD.af,
	       featureData=fD.af,
	       lowerDetectionLimit=1,
	       expressionFamily=negbinomial.size())

# Trajectory inference according to vignette
phenoData(cds)$Size_Factor <- pD$sf
cds <- estimateDispersions(cds)

# Feature Selection
genes <- rownames(hvg.out)
cds <- setOrderingFilter(cds,genes)

# Dim Reduction 
cds <- reduceDimension(cds, max_components=2,norm_method="log")
cds <- orderCells(cds,reverse=TRUE)

# Plot trajectory colored by clusters
p0 <- plot_cell_trajectory(cds,x=1,y=2, color_by="SuperCluster", cell_size=1,show_branch_points=FALSE) 

# Save
pD.monoc <- pData(cds)
monoc <- list("pD"=pD.monoc,
	      "plot"=p0)
saveRDS(monoc,"../data/Robjects/secondRun_2500/ExpressionList_Monocle.rds")

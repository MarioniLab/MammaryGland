# Differentation inference with monocle

library(plyr)
library(dplyr)
library(ggplot2)
library(monocle)
library(RColorBrewer)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# ---- Prepocessing ----

# Cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier & pD$Condition %in% c("NP","G") &
	!pD$cluster %in% c(6,7,9)
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep <- rowMeans(m)>0.1
m <- m[keep,]
fD <- fD[keep,]

# Normalize
m.norm <- t(t(m)/pD$sf)
rownames(m.norm) <- fD$symbol

# High Var
brennecke <- BrenneckeHVG(m.norm,fdr=0.1)
fD$highVar <- fD$id %in% brennecke

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
genes <- filter(fD, highVar) %>% .$id
cds <- setOrderingFilter(cds,genes)

# Dim Reduction 
cds <- reduceDimension(cds, max_components=2,norm_method="log")
cds <- orderCells(cds,reverse=TRUE)

# Plot trajectory colored by clusters
pal <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,8,9)]
p0 <- plot_cell_trajectory(cds,x=1,y=2, color_by="cluster", cell_size=1,show_branch_points=FALSE) +
    scale_color_manual(values=pal)

# Save
pD.monoc <- pData(cds)
monoc <- list("pD"=pD.monoc,
	      "plot"=p0)
saveRDS(monoc,"../data/Robjects/ExpressionList_Monocle.rds")

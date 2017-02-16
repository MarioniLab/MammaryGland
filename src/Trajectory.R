library(scran)
library(dplyr)
library(ggplot2)
library(monocle)
source("functions.R")

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier #& (pD$cluster %in% c(6,7,10)) & pD$Condition %in% c("V","P","L","I")
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep1 <- rowSums(m!=0)>10
keep2 <- rowSums(m)>20
keep <- keep1 & keep2
m <- m[keep,]
fD <- fD[keep,]

#Normalize
m.norm <- t(t(m)/pD$sf)

#High Var
brennecke <- BrenneckeHVG(m.norm,fdr=0.1)
fD$highVar <- fD$id %in% brennecke



pD.af <- new("AnnotatedDataFrame", data=pD)
sampleNames(pD.af) <- pD$barcode
fD.af <- new("AnnotatedDataFrame", data=fD)
sampleNames(fD.af) <- fD$id

#Initiate Monocle
cds <- newCellDataSet(as.matrix(m),
	       phenoData=pD.af,
	       featureData=fD.af,
	       lowerDetectionLimit=1,
	       expressionFamily=negbinomial.size())
phenoData(cds)$Size_Factor <- pD$sf
cds <- estimateDispersions(cds)

#Feature Selection
genes <- filter(fD, highVar) %>% .$id
cds <- setOrderingFilter(cds,genes)

#
cds <- reduceDimension(cds, max_components=2)
cds <- orderCells(cds,reverse=FALSE)

plot_cell_trajectory(cds,x=1,y=2, color_by="Pseudotime", show_branch_points=TRUE, cell_size=3) + facet_wrap(~cluster) + scale_color_brewer(palette="Paired") 

# gns <- filter(fD, symbol %in% c("Elf5","Prlr","Krt18","Krt5","Hp")) %>% .$id
# plot_genes_in_pseudotime(cds[gns,], color_by="cluster")
# plot_genes_branched_pseudotime(cds[gns,], color_by="cluster",branch_point=1) + scale_color_brewer(palette="Paired") 



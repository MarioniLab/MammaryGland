library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
source(file.path("functions.R"))

rnd_seed <- 300
dataList <- readRDS(file.path("../data/Robjects/ExpressionList_Clustered.rds"))
deList <- readRDS(file.path("../data/Robjects/DEList.rds"))
deMethGenes <- read.table("../data/methylationData/HypoMethInParous.txt", header=TRUE, stringsAsFactor=FALSE) %>% .$Gene


#Commence
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# Cells
keepCells <- pD$PassAll & ! (pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep1 <- rowSums(m!=0)>10
keep2 <- rowSums(m)>20
keep <- keep1 & keep2
m <- m[keep,]
fD <- fD[keep,]

#Normalization
m <- t(t(m)/pD$sf)

#FCs
genes <- c("B4galt1",
	   "Bax",
	   "Stat5a",
	   "Vdr",
	   "Atp7b",
	   "Med16",
	   "Socs2",
	   "Med13",
	   "Hif1a",
	   "Med12l",
	   "Med13l",
	   "Eif2ak3",
	   "Agpat6",
	   "Csn3",
	   "Med18",
	   "Cdo1",
	   "Xdh",
	   "Vegfa",
	   "Stat5b",
	   "Med1",
	   "Birc2")

genes <- filter(fD, symbol %in% genes) %>% .$id

deResults <- deList[[4]]
forPlot <- deResults[,c("Gene","logFC.vs.5","top.vs.5")]
forPlot <- mutate(forPlot, isHypoMeth= Gene %in% deMethGenes) %>% 
	   rename(foldChange=logFC.vs.5)

x <- filter(forPlot, foldChange >0, top.vs.5 < 1000)
y <- filter(forPlot, foldChange <0, top.vs.5 < 1000)

p <- ggplot(forPlot, aes(y=foldChange, x=isHypoMeth)) +
    geom_boxplot() +
    theme_bw()

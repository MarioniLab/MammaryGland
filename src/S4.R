# S5

library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(gridGraphics)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# add SubClusterNumbers
dataList2 <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
pD.add <- dataList2[[2]]
rm(dataList2)
pD.add <- pD.add[,c("barcode","SubClusterNumbers")]
pD <- left_join(pD,pD.add,by="barcode")

# Remove QC fail 
m <- m[,pD$PassAll]
pD <- pD[pD$PassAll,]

m <- t(t(m)/pD$sf)


# First step
p1 <- ggplot(pD,aes(tSNE1,tSNE2,color=Cluster)) +
    geom_point()

# Second round
fp2 <- pD
fp2$SubCluster <- as.character(fp2$SubCluster)
fp2$SubCluster <- substr(fp2$SubCluster,nchar(fp2$SubCluster),nchar(fp2$SubCluster))
p2 <- ggplot(fp2,aes(tSNE1,tSNE2,color=SubCluster)) +
    geom_point(size=.8) +
    facet_wrap(~Cluster)

subP1 <- plot_grid(p1,p2, labels=c("a","b"))

# ---- PlotImmuneCellMarkers ----

genes <- c("Cd52","Cd74","Cd72")
p3 <- plotGeneDist(m, pD, fD, genes, colorBy="SubClusterNumbers")
p3 <- p3 %+% facet_grid(variable~.) %+% xlab("Cluster") %+% ylab("Log-Expression") %+% 
    theme(axis.text.x=element_text(angle=45,hjust=1))


# ---- PlotFibroblastMarkers ----

genes <- c("Col3a1","Col5a1","Col6a1","Fn1")
p4 <- plotGeneDist(m, pD, fD, genes, colorBy="SubClusterNumbers")
p4 <- p4 %+% facet_grid(variable~.) %+% xlab("Cluster") %+% ylab("Log-Expression") %+% 
    theme(axis.text.x=element_text(angle=45,hjust=1))
# ---- PlotEnodthelialCellMarkers ----

genes <- c("Eng","S1pr1","Emcn")
p5 <- plotGeneDist(m, pD, fD, genes, colorBy="SubClusterNumbers")
p5 <- p5 %+% facet_grid(variable~.) %+% xlab("Cluster") %+% ylab("Log-Expression") %+% 
    theme(axis.text.x=element_text(angle=45,hjust=1))

# ---- PlotPericyteMarkers ----

genes <- c("Pdgfrb","Cspg4","Anpep","Des")
p6 <- plotGeneDist(m, pD, fD, genes, colorBy="SubClusterNumbers")
p6 <- p6 %+% facet_grid(variable~.) %+% xlab("Cluster") %+% ylab("Log-Expression") %+% 
    theme(axis.text.x=element_text(angle=45,hjust=1))


subP2 <- plot_grid(p3,p4,p5,p6, nrow=1, labels=c("c","d","e","f"))

fullP <- plot_grid(subP1,subP2,nrow=2)
cairo_pdf("../paper/figures/S4.pdf",height=12.41,width=17.54)
fullP
dev.off()

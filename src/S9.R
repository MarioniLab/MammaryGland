# S9

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
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# Remove outlier and immune cells
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
p3 <- plotGeneDist(m, pD, fD, genes, colorBy="SubCluster")
p3 <- p3 %+% facet_grid(variable~.) %+% xlab("SubCluster") %+% ylab("Log-Expression")

# ---- PlotEnodthelialCellMarkers ----

genes <- c("Eng","S1pr1","Emcn")
p4 <- plotGeneDist(m, pD, fD, genes, colorBy="SubCluster")
p4 <- p4 %+% facet_grid(variable~.) %+% xlab("Cluster") %+% ylab("Log-Expression")

# ---- PlotPericyteMarkers ----

genes <- c("Pdgfrb","Cspg4","Anpep","Des")
p5 <- plotGeneDist(m, pD, fD, genes, colorBy="SubCluster")
p5 <- p5 %+% facet_grid(variable~.) %+% xlab("Cluster") %+% ylab("Log-Expression")

subP2 <- plot_grid(p3,p4,p5, nrow=1, labels=c("c","d","e"))

fullP <- plot_grid(subP1,subP2,nrow=2)
cairo_pdf("../paper/figures/S9.pdf",height=12.41,width=17.54)
fullP
dev.off()

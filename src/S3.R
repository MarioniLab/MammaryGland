# Figure S1
library(scran)
library(plyr)
library(dplyr)
library(knitr)
library(ggplot2)
library(gridExtra)
library(Rtsne)
library(cowplot)
source("functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC.rds")
sumDat <- read.csv("../data/CellRangerData/secondRun_2500/Summary.csv")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]

# Summary table 
sumry <- group_by(pD, SampleID) %>%
    summarize("Number of cells"=n(),
	      "Total molecules"=median(UmiSums),
	      "Genes Detected"=median(GenesDetected))
sumry <- left_join(sumry, sumDat[,c("SampleID","NumberOfReads","Saturation")])

p1 <- tableGrob(sumry,rows=NULL,)

# Load plots from QCAnalysis
# Illustrate thresholds
gdHist <- ggplot(pD, aes(x=GenesDetected,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=100) +
    scale_x_log10() +
    xlab("Total number of genes detected") +
    facet_wrap(~Condition) 

libSizeHist <- ggplot(pD, aes(x=UmiSums,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=100) +
    scale_x_log10() +
    facet_wrap(~Condition) +
    xlab("Total number of unique molecules") 

cellViability <- ggplot(pD, aes(x=prcntMito, y=GenesDetected, color=Condition, shape=Replicate))+
    geom_point() +
    xlab("Percentage of Mitochondrial RNA molecules") +
    ylab("Total number of genes detected")

cairo_pdf("../paper/figures/S3.pdf",height=12.41,width=17.54)
plot_grid(p1,gdHist,libSizeHist,cellViability, labels="auto")
dev.off()

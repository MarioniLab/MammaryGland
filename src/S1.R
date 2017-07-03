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
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList.rds")
sumDat <- read.csv("../data/CellRangerData/secondRun_2500/Summary.csv")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]

pD$UmiSums<- colSums(m)
pD$GenesDetected <- colSums(m!=0)

# Summary table 
sumry <- group_by(pD, SampleID) %>%
    summarize("Number of cells"=n(),
	      "Total molecules"=median(UmiSums),
	      "Genes Detected"=median(GenesDetected))
sumry <- left_join(sumry, sumDat[,c("SampleID","NumberOfReads","Saturation")])

p1 <- tableGrob(sumry,rows=NULL,)

# Load plots from QCAnalysis
source("QCAnalysis.R")
# dev.off()
# cairo_pdf("../paper/figures/S1.pdf",height=12.41,width=17.54)
plot_grid(p1,gdHist,libSizeHist,cellViability, labels="auto")
# dev.off()

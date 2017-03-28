library(scran)
library(plyr)
library(dplyr)
library(knitr)
library(ggplot2)
library(gridExtra)
library(Rtsne)
library(cowplot)
source("functions.R")
dataList <- readRDS("../data/Robjects/ExpressionList.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
#relevel pD Condition to get it in logical sequence
pD$Condition <- mapvalues(pD$Condition, from=c("V","P","L","I"),
			  to=c("NP", "G",
			       "L", "PI"))
pD$Condition <- factor(pD$Condition, levels=c("NP","G","L","PI"))
pD$SampleID <- mapvalues(pD$SampleID, from=c("V1","V2","P1","P2","L1","L2","I1","I2"),
			  to=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2"))
pD$SampleID <- factor(pD$SampleID,levels=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2")) 

pD$UmiSums<- colSums(m)
pD$GenesDetected <- colSums(m!=0)

sumry <- group_by(pD, SampleID) %>%
    summarize("Number of cells"=n(),
	      "Total molecules"=median(UmiSums),
	      "Genes Detected"=median(GenesDetected)) %>%
mutate("Number of reads"=c(257639,124686,574902,372905,106789,111926,264971,914353)) %>%
mutate("Sequencing Saturation"=c(86.3,75.3,91.4,88.3,97.7,98,90.2,97.2)) %>%
rename(Sample=SampleID)

p1 <- tableGrob(sumry,rows=NULL,)


mMito <- m[fD$Mitochondrial,]
pD$prcntMito <- colSums(mMito)/colSums(m)
cellViability <- ggplot(pD, aes(x=prcntMito, y=GenesDetected, color=Condition))+
    geom_point() +
    theme_bw()
leftmad <- function(x) {
    m <- median(x)
    dev <- abs(x-m)
    leftmad <- median(dev[x<=m])
    return(leftmad)
}
smryByGroup <- group_by(pD,Condition) %>%
    summarize(medianPrcntMito=median(prcntMito),
	      madPrcntMito=mad(prcntMito),
	      threshold_PrcntMito=0.05,
	      mGenesDetected=median(log10(GenesDetected)),
	      madGenesDetected=leftmad(log10(GenesDetected)),
	      threshold_GenesDetected=max(mGenesDetected-4*madGenesDetected,log10(500)),
	      mUmiSums=median(log10(UmiSums)),
	      madUmiSums=leftmad(log10(UmiSums)),
	      threshold_UmiSums=max(mUmiSums-4*madUmiSums,log10(1000))) %>%
    select(Condition,starts_with("threshold")) 

#simple thresholds
MitoCutOff <- 0.05
#inititae Df
pD <- mutate(pD,
	     ThresholdViability = 0,
	     ThresholdGenesDet = 0,
	     ThresholdLibSize = 0)

grps <- as.character(unique(pD$Condition))
for (grp in grps) {
    thrs <- filter(smryByGroup, Condition==grp) %>% select(-Condition) %>% t() %>% as.vector()
    names(thrs) <- filter(smryByGroup, Condition==grp) %>% select(-Condition) %>% t() %>% rownames()
    pD <- mutate(pD,
		 ThresholdViability= ifelse(Condition==grp, MitoCutOff, ThresholdViability),
		 ThresholdGenesDet= ifelse(Condition==grp, 10^thrs["threshold_GenesDetected"],ThresholdGenesDet),
		 ThresholdLibSize= ifelse(Condition==grp, 10^thrs["threshold_UmiSums"],ThresholdLibSize))
}

pD <- mutate(pD,
	     PassViability=prcntMito < ThresholdViability,
	     PassGenesDet=GenesDetected > ThresholdGenesDet,
	     PassLibSize=UmiSums > ThresholdLibSize,
	     PassAll= PassViability & PassGenesDet & PassLibSize)

gdHist <- ggplot(pD, aes(x=GenesDetected,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=50) +
    geom_vline(data=smryByGroup,aes(xintercept=10^threshold_GenesDetected),color="red",lty="longdash") +
    scale_x_log10() +
    xlab("Number of Genes Detected") +
    facet_wrap(~Condition) 

libSizeHist <- ggplot(pD, aes(x=UmiSums,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=50) +
    geom_vline(data=smryByGroup,aes(xintercept=10^threshold_UmiSums),color="red",lty="longdash") +
    scale_x_log10() +
    xlab("Total Number of Molecules") +
    facet_wrap(~Condition) 

cellViability <- cellViability %+% pD
cellViability <- cellViability + aes(color=PassViability) +
    annotate("rect",ymin=-Inf, ymax=Inf, xmax=Inf, xmin=MitoCutOff,
	     fill="grey", alpha=0.3) +
scale_color_manual(values=c("red","black")) +
xlab("Percent of Molecules from Mitochondrial transcripts")+
ylab("Number of Genes Detected")

cairo_pdf("../paper/figures/S1.pdf",height=12.41,width=17.54)
plot_grid(p1,gdHist,libSizeHist,cellViability, labels="auto")
dev.off()

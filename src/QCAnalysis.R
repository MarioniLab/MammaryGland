####################
#
#QC-Analysis
#
###################

library(scran)
library(dplyr)
library(knitr)
library(ggplot2)
library(Rtsne)
library(cowplot)
source("functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/ExpressionList.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rnd_seed <- 300

# ---- QCOverview ----

# Number of Cells
table(pD$SampleID)
table(pD$Condition)

# Sequencing Depth and Genes detected
pD$UmiSums<- colSums(m)
pD$GenesDetected <- colSums(m!=0)
genesDetected <- ggplot(pD, aes(x=SampleID,y=GenesDetected,fill=Condition)) +
    geom_violin(draw_quantiles=0.5)+
    scale_y_log10() +
    ylab("Total number of genes detected") +
    theme_bw()
LibrarySize <- ggplot(pD, aes(x=SampleID,y=UmiSums,fill=Condition)) +
    geom_violin(draw_quantiles=0.5)+
    scale_y_log10() +
    ylab("Total number of molecules") +
    theme_bw()

# Cell Viability
mMito <- m[fD$Mitochondrial,]
pD$prcntMito <- colSums(mMito)/colSums(m)
cellViability <- ggplot(pD, aes(x=prcntMito, y=GenesDetected, color=Condition, shape=Replicate))+
    geom_point() +
    theme_bw()

# Genewise QC
ntop <- 50 
mRel <- t(t(m)/colSums(m))
rownames(mRel)  <- fD$symbol
topExpressed <- rowMedians(mRel)
names(topExpressed) <- rownames(mRel)
topExpressed <- topExpressed %>% sort(.,decreasing=TRUE) %>% names
plotData <- t(mRel)[,topExpressed[1:ntop]] %>%
    reshape2::melt() %>%
    rename(Cell=Var1,
	   Gene=Var2,
	   RelativeExpression=value)
topGenes <- ggplot(plotData, aes(x=Gene, y=RelativeExpression)) +
    geom_boxplot() +
    coord_flip() +
    theme_bw()
freqOfExp <- m!=0
rownames(freqOfExp) <- fD$symbol
freqOfExp <- sort(rowSums(freqOfExp)/ncol(freqOfExp),decreasing=TRUE)
plotData <- data.frame("Gene"=names(freqOfExp),"Frequency"=freqOfExp)
topFreq <- ggplot(plotData[1:ntop,], aes(x=factor(Gene,levels=Gene),y=Frequency)) +
    geom_bar(stat="identity") +
    coord_flip() +
    xlab("Gene") +
    theme_bw()


plot_grid(genesDetected,LibrarySize,cellViability,
	  topGenes,topFreq)


# ---- QCThresholding ----


# Left MAD function for thresholding
leftmad <- function(x) {
    m <- median(x)
    dev <- abs(x-m)
    leftmad <- median(dev[x<=m])
    return(leftmad)
}

# Define hard thresholds
MitoCutOff <- 0.05
genesDetectedCutOff <- 500 # only applies to L sample
librarySizeCutOff <- 1000 # only applies to L sample

smryByGroup <- group_by(pD,Condition) %>%
    summarize(threshold_PrcntMito=0.05,
	      mGenesDetected=median(log10(GenesDetected)),
	      madGenesDetected=leftmad(log10(GenesDetected)),
	      threshold_GenesDetected=max(mGenesDetected-4*madGenesDetected,log10(genesDetectedCutOff)),
	      mUmiSums=median(log10(UmiSums)),
	      madUmiSums=leftmad(log10(UmiSums)),
	      threshold_UmiSums=max(mUmiSums-4*madUmiSums,log10(librarySizeCutOff))) %>%
    select(Condition,starts_with("threshold")) 
print(smryByGroup)


# initiate empty df
pD <- mutate(pD,
	     ThresholdViability = 0,
	     ThresholdGenesDet = 0,
	     ThresholdLibSize = 0)

# Thresholding per group
grps <- as.character(unique(pD$Condition))
for (grp in grps) {
    thrs <- filter(smryByGroup, Condition==grp) %>% select(-Condition) %>% t() %>% as.vector()
    names(thrs) <- filter(smryByGroup, Condition==grp) %>% select(-Condition) %>% t() %>% rownames()
    pD <- mutate(pD,
		 ThresholdViability= ifelse(Condition==grp, thrs["threshold_PrcntMito"], ThresholdViability),
		 ThresholdGenesDet= ifelse(Condition==grp, 10^thrs["threshold_GenesDetected"],ThresholdGenesDet),
		 ThresholdLibSize= ifelse(Condition==grp, 10^thrs["threshold_UmiSums"],ThresholdLibSize))
}

pD <- mutate(pD,
	     PassViability=prcntMito < ThresholdViability,
	     PassGenesDet=GenesDetected > ThresholdGenesDet,
	     PassLibSize=UmiSums > ThresholdLibSize,
	     PassAll= PassViability & PassGenesDet & PassLibSize)

gdHist <- ggplot(pD, aes(x=GenesDetected,y=..density..)) +
    geom_histogram(fill="white",color="black") +
    geom_vline(data=smryByGroup,aes(xintercept=10^threshold_GenesDetected),color="red",lty="longdash") +
    scale_x_log10() +
    xlab("Total number of genes detected") +
    facet_wrap(~Condition) +
    theme_bw()

libSizeHist <- ggplot(pD, aes(x=UmiSums,y=..density..)) +
    geom_histogram(fill="white",color="black") +
    geom_vline(data=smryByGroup,aes(xintercept=10^threshold_UmiSums),color="red",lty="longdash") +
    scale_x_log10() +
    facet_wrap(~Condition) +
    xlab("Total number of unique molecules") +
    theme_bw()

cellViability <- cellViability %+% pD
cellViability <- cellViability + aes(color=PassViability) +
    annotate("rect",ymin=-Inf, ymax=Inf, xmax=Inf, xmin=MitoCutOff,
	     fill="grey", alpha=0.3) 

# Overview over cells removed
table(pD$Condition,pD$PassGenesDet)
table(pD$Condition,pD$PassLibSize)
table(pD$Condition,pD$PassViability)
table(pD$Condition,pD$PassAll)

plot_grid(gdHist,libSizeHist,cellViability)

# Apply Thresholds
pD.filtered <- filter(pD, PassAll)
m.filtered <- m[,as.character(pD.filtered$barcode)]

# Thresholding on lowly expressed genes
isexpThreshold <- 10
expThreshold <- 50*isexpThreshold
keep1 <- rowSums(m.filtered!=0) > isexpThreshold
keep2 <- rowSums(m.filtered) > expThreshold
fD$keep <- keep1 & keep2
m.filtered <- m.filtered[fD$keep,]
fD.filtered <- fD[fD$keep,]

# Final matrix
dim(m.filtered)

# ---- NormAndHVG ----

# Estimate size factors using scran
clusters <- quickCluster(m.filtered)
pD.filtered$sf <- computeSumFactors(m.filtered,clusters=clusters)

plot(log10(colSums(m.filtered))~log10(pD.filtered$sf),main="Library Size versus Size Factors (Log10-Scale)",
     pch=20,xlab="Size Factors",ylab="Library Size")

# Normalize count matrix
m.norm <- t(t(m.filtered)/pD.filtered$sf)


# Highly variable genes
brennecke <- BrenneckeHVG(m.norm,fdr=0.1,minBiolDisp=0.25)
fD.filtered$highVar <- fD.filtered$id %in% brennecke


# ---- SaveData ----
pD.add <- select(pD.filtered, barcode, sf)
fD.add <- select(fD.filtered, id, highVar)
pD <- left_join(pD,pD.add, by="barcode")
fD <- left_join(fD,fD.add, by="id")
stopifnot(identical(as.character(pD$barcode),colnames(m)))
stopifnot(identical(as.character(fD$id),rownames(m)))
out <- list()
out[[1]] <- m
out[[2]] <- pD
out[[3]] <- fD
saveRDS(out,file="../data/Robjects/ExpressionList_QC.rds")

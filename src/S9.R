# Estimate size factors using scran
library(dplyr)
library(scran)
library(Rtsne)
library(knitr)
library(gridExtra)
source("functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
# set.seed(300)
# dataList <- subSample(dataList, cell.number=2000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
m <- t(t(m)/pD$sf)

# Doublet definition: below 7% in all samples where it was present
tab <- ftable(pD$SubClusterNumbers, pD$SampleID)
tab <- as.matrix(t(t(tab)/colSums(tab)))
pot.dob <- rownames(tab)[which(rowMax(tab)<0.07)]

freqDf <- data.frame(as.table(tab))
ordr <- group_by(freqDf, Var1) %>%
    summarize(tot=sum(Freq)) %>%
    dplyr::arrange(desc(tot)) %>%
    .$Var1 %>%
    as.character()

freqDf$Var1 <- factor(freqDf$Var1,levels=ordr)
colnames(freqDf) <- c("Cluster","Sample","Frequency")

p1 <- ggplot(freqDf, aes(x=Cluster,y=Frequency, fill=Sample)) +
    geom_bar(stat="identity",position="dodge") +
    geom_hline(yintercept=0.07,color="red",lty="dashed") +
    theme(axis.text.x=element_text(angle=45, vjust=0.5, size=16))

# for Table

out <- data.frame()
smpls <- levels(pD$SampleID)
for (smpl in smpls) {
    pD.sub <- pD[pD$SampleID==smpl,]
    keepC <- names(table(pD.sub$SubClusterNumbers)[which(table(pD.sub$SubClusterNumbers)>1)])
    pD.sub <- pD.sub[pD.sub$SubClusterNumbers %in% keepC,]
    m.sub <- m[,pD.sub$barcode]
    cluster <- factor(pD.sub$SubClusterNumbers)
    tmp <- data.frame(numeric(nrow(m)),check.names=FALSE)
    colnames(tmp) <- cluster[1]
    for (clust in levels(cluster)) { 
	expr <- rowMeans(as.matrix(log2(m.sub[,pD.sub$SubClusterNumbers==clust]+1)))
	colname <- clust
	tmp[,colname] <- expr
    }
    tmp <- cor(tmp) # compute correlation
    tmp <- apply(tmp,1, function(x) sum(x>0.9)) #how many are more than 0.9
    tmp <- names(tmp[which(tmp>=3)])
    tmp1 <- data.frame("potDoublet"=pot.dob,
		       "Sample"=smpl)
    tmp1$presentInSample <- tmp1$potDoublet %in% keepC
    tmp1$CorWithAtLeastTwo <- tmp1$potDoublet %in% tmp
    out <- rbind(out,tmp1)
}
out <- out[out$presentInSample,]

xout <- data.frame()
for (smpl in c("G1","G2")) {
    pD.sub <- pD[pD$SampleID==smpl,]
    m.sub <- m[,pD.sub$barcode]
    cluster <- factor(pD.sub$SubClusterNumbers)
    tmp <- data.frame(numeric(nrow(m)),check.names=FALSE)
    colnames(tmp) <- cluster[1]
    for (clust in levels(cluster)) { 
	expr <- rowMeans(as.matrix(log2(m.sub[,pD.sub$SubClusterNumbers==clust]+1)))
	colname <- clust
	tmp[,colname] <- expr
    }
    x1 <- data.frame("DoubletCluster"=tmp[,"C20"],
		     "Singleton"=tmp[,"C8"],
		     "Cluster"="C8")
    x2 <- data.frame("DoubletCluster"=tmp[,"C20"],
		     "Singleton"=tmp[,"C13"],
		     "Cluster"="C13")
    x2 <- rbind(x1,x2)
    x3 <- data.frame("DoubletCluster"=tmp[,"C20"],
		     "Singleton"=rowMeans(tmp[,c("C13","C8")]),
		     "Cluster"="Mean(C8,C13)")
    x4 <- rbind(x2,x3)
    x4$Sample <- paste0("Sample-",smpl)
    xout <- rbind(xout,x4)
}

p2 <- ggplot(xout, aes(Singleton, DoubletCluster)) +
    geom_point() +
    geom_abline(slope=1, intercept=0, color="red", lty="dashed") +
    facet_grid(Sample~Cluster) +
    xlab("Log-Expression") +
    ylab("Log-Expression") +
    theme_bw() +
    ylab("Log-Expression of C20")

pD <- group_by(pD, SampleID) %>%
    dplyr::mutate(GenesDetected=GenesDetected/median(GenesDetected)) %>%
    dplyr::mutate(UmiSums=UmiSums/median(UmiSums)) %>%
    ungroup()

ordr <- group_by(pD, SubClusterNumbers) %>%
    summarize(tot=median(GenesDetected))  %>%
    dplyr::arrange(desc(tot)) %>%
    .$SubClusterNumbers %>%
    as.character()

fp3 <- pD
fp3$SubClusterNumbers <- factor(fp3$SubClusterNumbers, levels=ordr)

p3 <- ggplot(fp3, aes(y=GenesDetected, x=SubClusterNumbers)) +
    geom_boxplot() +
    scale_y_log10() +
    theme_bw() +
    ylab("Total number of genes detected") +
    xlab("Cluster")

fp4 <- pD

ordr <- group_by(pD, SubClusterNumbers) %>%
    summarize(tot=median(UmiSums))  %>%
    dplyr::arrange(desc(tot)) %>%
    .$SubClusterNumbers %>%
    as.character()

fp4$SubClusterNumbers <- factor(fp4$SubClusterNumbers, levels=ordr)

p4 <- ggplot(fp4, aes(y=UmiSums, x=SubClusterNumbers)) +
    geom_boxplot() +
    scale_y_log10() +
    theme_bw() +
    ylab("Total number of molecules detected") +
    xlab("Cluster")

library(cowplot)


cairo_pdf("../paper/figures/S9.pdf",height=12.41,width=17.54)
plot_grid(p1,p2,p3,p4, labels="auto")
dev.off()

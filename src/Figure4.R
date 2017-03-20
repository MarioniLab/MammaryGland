library(scran)
library(plyr)
library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
library(dynamicTreeCut)
library(Rtsne)
library(reshape2)
library(pheatmap)
source("functions.R")

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- pD$PassAll & !(pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep <- rowMeans(m) > 0.1
m <- m[keep,]
fD <- fD[keep,]

#Normalize
m.norm <- t(t(m)/pD$sf)

#HVG
brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)

#Rename condition and cluster
pD$Condition <- mapvalues(pD$Condition, from=c("V","P","L","I"),
			  to=c("Virgin", "Pregnancy",
			       "Lactation", "Post-Involution"))

pD$cluster <- mapvalues(pD$cluster, from=c(1,2,3,4,5,6,7,9,10),
			  to=c(1,2,3,4,5,6,7,8,9))


out <- data.frame("Cluster 1"=numeric(nrow(m)))

for (clust in levels(pD$cluster)[-1]) {
    expr <- rowMeans(log2(m.norm[,pD$cluster==clust]+1))
    colname <- paste0("Cluster.",clust)
    out[,colname] <- expr
}
library(pheatmap)

meanSim <- cor(out,method="spearman")
simMat <- pheatmap(meanSim)
dev.off()

######################
#
#DE-Analysis
#
######################

m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- pD$PassAll & !(pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

#Rename condition and cluster
pD$Condition <- mapvalues(pD$Condition, from=c("V","P","L","I"),
			  to=c("Virgin", "Pregnancy",
			       "Lactation", "Post-Involution"))

pD$cluster <- mapvalues(pD$cluster, from=c(1,2,3,4,5,6,7,9,10),
			  to=c(1,2,3,4,5,6,7,8,9))


comps <- c("Virgin","Post-Involution")
out <- list()
for (comp in comps) {
    pD.sub <- filter(pD,Condition %in% comp, !(cluster %in% c(6,7,9)))
    m.sub <- m[,as.character(pD.sub$barcode)]
    keep <- rowMeans(m.sub) > 0.01
    m.sub <- m.sub[keep,]
    fD.sub <- fD[keep,]
    rownames(m.sub) <- fD.sub$symbol

    library(edgeR)
    nf <- log(pD.sub$sf/pD.sub$UmiSums)
    pD.sub$nf <- exp(nf-mean(nf))
    y <- DGEList(counts=m.sub,
		 samples=pD.sub,
		 genes=fD.sub,
		 norm.factors=pD.sub$nf)


    choice <- ifelse(comp=="Virgin",5,4)
    cluster <- factor(as.numeric(pD.sub$cluster==choice))
    de.design <- model.matrix(~cluster)
    y <- estimateDisp(y, de.design, prior.df=0,trend="none")
    fit <- glmFit(y, de.design)
    res <- glmTreat(fit)
    xxx <- topTags(res,n=Inf,sort.by="PValue")
    topTab <- xxx$table
    out[[comp]] <- topTab
}
tabPar <- out[["Virgin"]][1:500,]
tabNulPar <- filter(out[["Post-Involution"]],symbol %in% tabPar$symbol)
rownames(tabPar) <- tabPar$symbol
rownames(tabNulPar) <- tabNulPar$symbol
tabNulPar <- tabNulPar[tabPar$symbol,]

forPlot <- data.frame("NullParFC"=tabPar$logFC,
		      "ParousFC"=tabNulPar$logFC,
		      "Gene"=tabPar$symbol)

library(cowplot)

interest <- filter(forPlot, Gene %in% c("Aldh1a3","Elf5","Hey1","Sox10","Cd14",
					"Prlr","Esr1","Pgr"))
p <- ggplot(forPlot, aes(x=NullParFC,y=ParousFC)) +
    geom_point(color="grey50") +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), size=4, color="red",pch=1) +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), color="black") +
    geom_label(data=interest, aes(x=NullParFC,y=ParousFC,label=Gene),nudge_x=0.7) +
    geom_hline(yintercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    coord_equal(xlim=c(-5,5),ylim=c(-5,5))



subP0 <- plot_grid(simMat[[4]],p,labels=c("a","b"))


######
#
#
#Hypo-Methylated Genes
#
#
#######
deMethGenes <- read.table("../data/methylationData/HypoMethInParous.txt", header=TRUE, stringsAsFactor=FALSE) %>% .$Gene
deList <- readRDS(file.path("../data/Robjects/DEList.rds"))


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
deRes <- deResults[,c("Gene","logFC.vs.5","top.vs.5")]
deRes <- mutate(deRes, isHypoMeth= Gene %in% deMethGenes) %>% 
	   rename(foldChange=logFC.vs.5)


p <- ggplot(filter(deRes, top.vs.5 <=1000), aes(y=foldChange, x=isHypoMeth)) +
    geom_violin(draw_quantiles=0.5) +
    xlab("Hypo-Methylated") +
    ylab("Log Fold Change") 


comps <- c("Up-regulated vs. 5", "Down-regulated vs. 5")
out <- list()
for (comp in comps) {
    if (comp=="Up-regulated vs. 5") {
    deDat <- filter(deRes, foldChange >0)
    } else{
    deDat <- filter(deRes, foldChange <0)
    }
    ranks <- rank(deDat$top.vs.5)
    ind <- deDat$isHypoMeth
    geneset <- ranks[ind]
    background <- ranks[!ind]
    ks <- ks.test(geneset,background,alternative="greater")

    forPlot <- data.frame("Ranks"=c(geneset,background),
			  "Set"=c(rep("Hypo-methylated",length(geneset)),
				  rep("Background",length(background))))
    
    out[[comp]] <- ggplot(forPlot, aes(x=Ranks,color=Set)) +
	stat_ecdf(geom="step",size=1.2) +
	annotate("text",x=2000,y=0.5,label=paste0("P=",round(ks$p.value,5))) +
	ggtitle(comp) +
	theme(legend.directio="horizontal",
	      legend.title=element_blank())
}
leg <- get_legend(out[[1]])
out[[1]] <- out[[1]] %+% guides(colour=FALSE) 
out[[2]] <- out[[2]] %+% guides(colour=FALSE) 
subP1a <- plot_grid(p,out[[1]],out[[2]],nrow=1,align="h")
subP1 <- plot_grid(subP1a,leg,rel_heights=c(1,0.1),nrow=2)
fullP <- plot_grid(subP0,subP1,labels=c("","c"),nrow=2)

cairo_pdf("Figure4.pdf",width=12.41,height=14.54)
fullP
dev.off()

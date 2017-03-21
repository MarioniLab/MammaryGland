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
library(RColorBrewer)

meanSim <- cor(out,method="spearman")
simMat <- pheatmap(meanSim,
		   color=colorRampPalette((brewer.pal(n=7,
							 name="Greys")))(100),
		   treeheight_row=0)
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

interest <- filter(forPlot, Gene %in% c("Aldh1a3","Elf5","Sox10","Cd14",
					"Prlr","Esr1","Pgr"))
p <- ggplot(forPlot, aes(x=NullParFC,y=ParousFC)) +
    geom_point(color="grey50") +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), size=4, color="red",pch=1) +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), color="black") +
    geom_label(data=filter(interest,NullParFC < 0), aes(x=NullParFC,y=ParousFC,label=Gene), nudge_x=-1) +
    geom_label(data=filter(interest,NullParFC > 0), aes(x=NullParFC,y=ParousFC,label=Gene), nudge_x=1.2) +
    xlab("Log(FC) against luminal cells in NP gland") +
    ylab("Log(FC) against luminal cells in P gland") +
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
			  "Geneset"=c(rep("Hypo-methylated",length(geneset)),
				  rep("Background",length(background))))
    
    out[[comp]] <- ggplot(forPlot, aes(x=Ranks,color=Geneset)) +
	stat_ecdf(geom="step",size=1.2) +
	annotate("text",x=2000,y=0.5,label=paste0("P=",round(ks$p.value,5))) +
	ggtitle(comp) +
	scale_colour_manual(values=c("grey50","black")) +
	ylab("Proportion of Genes") +
	xlab("P-Value Ranks") +
	theme(legend.directio="horizontal",
	      legend.title=element_blank())
}
leg <- get_legend(out[[1]])
out[[1]] <- out[[1]] %+% guides(colour=FALSE) 
out[[2]] <- out[[2]] %+% guides(colour=FALSE) 
subP1a <- plot_grid(p,out[[1]],out[[2]],nrow=1,align="h")
subP1 <- plot_grid(subP1a,leg,rel_heights=c(1,0.1),nrow=2)

####
#
#
#DE-Analysis 4 vs 5 
#
####

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

stopifnot(identical(rownames(m.norm),as.character(fD$id)))
rownames(m.norm) <- as.character(fD$symbol)
Milkgenes <- c("Csn2","Csn3","Lalba","Csn1s1","Csn1s2a","Mfge8")
ImmuneGenes <- c("Lcn2","B2m","Hk2","H2-K1","Ltf","Nfkb1")
geneL <- list(Milkgenes,ImmuneGenes)
for (i in seq_along(geneL)) {
    genes <- geneL[[i]]
    forPl <- data.frame(t(m.norm)[,c(genes)]+1,
			barcode=colnames(m.norm))
    add <- select(pD,barcode, cluster, Condition) %>%
	mutate(barcode=as.character(barcode))
    forPl <- left_join(forPl,add,by="barcode") %>%
	filter(cluster %in% c(5,4))

    forPl <- melt(forPl,id=c("barcode","cluster","Condition"))

    library(RColorBrewer)
    pal <- brewer.pal(n=9,name="Paired")[c(4,5)]
    ExpPlot <- ggplot(forPl, aes(y=value,x=cluster,color=cluster)) +
	geom_jitter(size=0.9) +
	stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
		     geom = "crossbar", width = 1,color="black") +
	facet_grid(~variable) +
	ylab("Log Expression") +
	xlab("") +
	theme(axis.line.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank(),
	      strip.background=element_blank(),
	      strip.text=element_text(face="bold"),
	      legend.position="bottom",
	      legend.direction="horizontal",
	      ) +
	scale_colour_manual(values=pal)+
	guides(colour = guide_legend(override.aes = list(size=3))) +
	scale_y_log10()

    out[[i]] <- ExpPlot
}

out[[1]] <- out[[1]] %+% guides(color=FALSE)
subP05 <- plot_grid(plotlist=out,nrow=2)
fullP <- plot_grid(subP0,subP05,subP1,labels=c("","c","d"),nrow=3)

cairo_pdf("Figure4.pdf",width=12.41,height=14.54)
fullP
dev.off()

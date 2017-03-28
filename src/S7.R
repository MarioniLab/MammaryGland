library(scran)
library(plyr)
library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
library(ggrepel)
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


out <- data.frame("C1"=numeric(nrow(m)))

for (clust in levels(pD$cluster)[-1]) {
    expr <- rowMeans(log2(m.norm[,pD$cluster==clust]+1))
    colname <- paste0("C",clust)
    out[,colname] <- expr
}
library(pheatmap)
library(RColorBrewer)

meanSim <- asSim(as.matrix(dist(t(out))))
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
tabNulPar <- out[["Virgin"]][1:500,]
tabPar <- filter(out[["Post-Involution"]],symbol %in% tabNulPar$symbol)
rownames(tabNulPar) <- tabNulPar$symbol
rownames(tabPar) <- tabPar$symbol
tabNulPar <- tabNulPar[tabPar$symbol,]

forPlot <- data.frame("NullParFC"=tabNulPar$logFC,
		      "ParousFC"=tabPar$logFC,
		      "Gene"=tabPar$symbol)

library(cowplot)
library(RColorBrewer)

#Try PCA instead
m.sub <- t(t(m.sub)/pD.sub$sf)
m.sub <- m.sub[tabNulPar$symbol,]
exps <- t(log2(m.sub+1))
pc <- prcomp(exps)
pD.sub$PC1 <- pc$x[,1]
pD.sub$PC2 <- pc$x[,2]
pal <- brewer.pal(name="Paired",n=9)[c(1,3,4,5,8)]

pcplot <- ggplot(pD.sub, aes(PC1,PC2,color=cluster)) +
    geom_point() +
    scale_color_manual(values=pal)

loadngs <- pc$rotation[,1]
loddf <- data.frame("symbol"=names(loadngs),
		    "Loadings"=loadngs)
ld <- join(tabNulPar,loddf) %>%
    mutate("Genes"=ifelse(logFC >0,"High in Progenitor","High in Differentiated"))
pLoad <- ggplot(ld, aes(x=Genes,y=Loadings)) +
    geom_boxplot() +
    ylab("PC1 Loadings")

subp0 <- plot_grid(pcplot, pLoad, rel_widths=c(1,.8),labels="auto") 



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

#########################################
    pD.sub <- filter(pD, (cluster %in% c(4,5)))
    m.sub <- m[,as.character(pD.sub$barcode)]
    keep <- rowMeans(m.sub) > 0.1
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


    choice <- 4
    cluster <- factor(as.numeric(pD.sub$cluster==choice))
    de.design <- model.matrix(~cluster)
    y <- estimateDisp(y, de.design, prior.df=0,trend="none")
    fit <- glmFit(y, de.design)
    res <- glmTreat(fit,lfc=1)
    xxx <- topTags(res,n=Inf,sort.by="PValue")
    topTab <- xxx$table

deGenes <- filter(topTab, FDR < 0.01, logFC > 0) %>%. $symbol
allGenes <- topTab$symbol

# GO Analysis
 require("topGO")
 library(org.Mm.eg.db)

# ---- Data ----
univrs <- allGenes
alG <- factor(as.numeric(univrs %in% deGenes))
names(alG) <- univrs

# ---- GOanalysis ----

## prepare Data for topGO
GO.data <- new("topGOdata", description="Lib GO",ontology="BP", allGenes=alG, 
	       annot=annFUN.org, mapping="org.Mm.eg.db",
	       nodeSize=5, ID="symbol")
result.classic <- runTest(GO.data, statistic="fisher")
output <- GenTable(GO.data, Fisher.classic=result.classic, orderBy="topgoFisher", topNodes=20, numChar=10000)

output$Term <- factor(output$Term, levels=rev(output$Term))
p1 <- ggplot(output, aes(x=Term, y=-log10(as.numeric(Fisher.classic)))) +
    geom_bar(stat="identity",color="black",fill="white") +
    coord_flip() +
    ylab("-log10(P)") +
    geom_hline(yintercept=3,lty="dashed") +
    xlab("GO-Term [BP]") 

subp1 <- plot_grid(subp0,p1,nrow=2,labels=c("","c"))
caseinPlot <- readRDS("../data/Robjects/CaseinsPlot.rds")
cairo_pdf("../paper/figures/S7.pdf",width=12.41,height=17.54)
plot_grid(subp1,caseinPlot,nrow=2,labels=c("","d"))
dev.off()

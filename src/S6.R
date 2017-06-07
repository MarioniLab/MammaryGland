# S6

library(plyr)
library(dplyr)
library(ggplot2)
library(Rtsne)
library(reshape2)
library(edgeR)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
source("functions.R")

# Load Data
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

# This DEAnalysis differs to the ProgenitorDE.R script, progenitors are compared to luminals from a specific timepoint
comps <- c("NP","PI")
out <- list()
for (comp in comps) {
    pD.sub <- filter(pD,Condition %in% comp, !(cluster %in% c(6,7,9)))
    m.sub <- m[,as.character(pD.sub$barcode)]
    keep <- rowMeans(m.sub) > 0.1
    m.sub <- m.sub[keep,]
    fD.sub <- fD[keep,]
    rownames(m.sub) <- fD.sub$symbol

    nf <- log(pD.sub$sf/pD.sub$UmiSums)
    pD.sub$nf <- exp(nf-mean(nf))
    y <- DGEList(counts=m.sub,
		 samples=pD.sub,
		 genes=fD.sub,
		 norm.factors=pD.sub$nf)


    choice <- ifelse(comp=="NP",5,4)
    cluster <- factor(as.numeric(pD.sub$cluster==choice))
    de.design <- model.matrix(~cluster)
    y <- estimateDisp(y, de.design, prior.df=0,trend="none")
    fit <- glmFit(y, de.design)
    res <- glmTreat(fit)
    resTab <- topTags(res,n=Inf,sort.by="PValue")
    topTab <- resTab$table
    out[[comp]] <- topTab
}
tabNulPar <- out[["NP"]][1:500,]

# ---- PCA ----
genes <- tabNulPar$symbol %in% rownames(m.sub)
m.sub <- t(t(m.sub)/pD.sub$sf)
m.sub <- m.sub[genes,]
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


# ---- DEAnalysis4vs5 ----
topTab <- read.csv("../data/Robjects/C4vsC5DE.csv")

deGenes <- filter(topTab, FDR < 0.01, logFC > 0) %>%. $symbol
allGenes <- topTab$symbol

# ---- Data ----
univrs <- allGenes
alG <- factor(as.numeric(univrs %in% deGenes))
names(alG) <- univrs

# ---- GOanalysis ----
library("topGO")
library(org.Mm.eg.db)

# prepare Data for topGO
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
# cairo_pdf("../paper/figures/S6.pdf",width=9.92,height=14.028)
plot_grid(subp1,nrow=2,labels=c("","d"))
# dev.off()

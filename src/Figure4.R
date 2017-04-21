# Figure 4
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
library(RColorBrewer)
library(edgeR)
library(cowplot)
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

pD$cluster <- factor(pD$cluster) #drop unused levels
out <- data.frame("C1"=numeric(nrow(m)))
for (clust in levels(pD$cluster)) {
    expr <- rowMeans(log2(m.norm[,pD$cluster==clust]+1))
    colname <- paste0("C",clust)
    out[,colname] <- expr
}

meanSim <- asSim(as.matrix(dist(t(out))))
simMat <- pheatmap(meanSim,
		   color=colorRampPalette((brewer.pal(n=7,
							 name="Greys")))(100),
		   treeheight_row=0)
dev.off()


# ---- DEProgenitorVsRest ----
#reload data
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

keepCells <- pD$PassAll & !(pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

comps <- c(5,4)
out <- list()
for (choice in comps) {
    rmClust <- setdiff(c(5,4),choice)
    pD.sub <- filter(pD, !(cluster %in% c(6,7,9,rmClust)))
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


    cluster <- factor(as.numeric(pD.sub$cluster==choice))
    de.design <- model.matrix(~cluster)
    y <- estimateDisp(y, de.design, prior.df=0,trend="none")
    fit <- glmFit(y, de.design)
    res <- glmTreat(fit,lfc=1)
    resTab <- topTags(res,n=Inf,sort.by="PValue")
    topTab <- resTab$table
    out[[paste0("C",choice)]] <- topTab
}
tabNulPar <- out[["C5"]][1:500,]
tabPar <- filter(out[["C4"]],symbol %in% tabNulPar$symbol)
rownames(tabNulPar) <- tabNulPar$symbol
rownames(tabPar) <- tabPar$symbol
tabNulPar <- tabNulPar[tabPar$symbol,]

forPlot <- data.frame("NullParFC"=tabNulPar$logFC,
		      "ParousFC"=tabPar$logFC,
		      "Gene"=tabPar$symbol)


#genes to highlight
interest <- filter(forPlot, Gene %in% c("Kit","Hey1","Cd14",
					"Prlr","Esr1","Pgr"))

#FC plot
p <- ggplot(forPlot, aes(x=NullParFC,y=ParousFC)) +
    geom_point(color="grey50",size=2) +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), color="black") +
    geom_label_repel(data=interest, aes(x=NullParFC,y=ParousFC,label=Gene)) +
    xlab("LFC of C5 vs. luminal cells") +
    ylab("LFC of C4 vs. luminal cells") +
    geom_hline(yintercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    coord_equal(xlim=c(-5,5),ylim=c(-5,5))


#Combine a and b
subP0 <- plot_grid(simMat[[4]],p,labels=c("a","b"))

# ---- DEC4vsC5 ----

# reload data
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

keepCells <- pD$PassAll & !(pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

#prepare data for DE
pD.sub <- filter(pD, (cluster %in% c(4,5)))
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


choice <- 4
cluster <- factor(as.numeric(pD.sub$cluster==choice))
de.design <- model.matrix(~cluster)
y <- estimateDisp(y, de.design, prior.df=0,trend="none")
fit <- glmFit(y, de.design)
res <- glmTreat(fit,lfc=1)
resTab <- topTags(res,n=Inf,sort.by="PValue")
topTab <- resTab$table

#Highlight genes with lactation/immune annotation
lac <- c("Btn1a1","Lalba","B4galt1","Csn3","Csn1s2a","Csn1s1","Csn2","Hk2","Xdh","Vegfa")
immuno <- c("Hp","Slpi","H2-K1", "B2m", "H2-Q7", "Lbp", "Tlr2", "Ltf", "Ifit1",
	    "Cd1d1")

#Write DE table for supps
forxls <- select(topTab, id, symbol, logFC, unshrunk.logFC, logCPM, PValue, FDR)
# write.csv(forxls,file="../paper/supps/DE_C4vsC5.csv",quote=FALSE)

#Highlight top DE genes in Volcano plot
topUp <- filter(topTab, FDR < 0.01) %>%
    arrange(logFC) %>% .$symbol %>% as.character() %>% .[1:5]
topDown <- filter(topTab, FDR < 0.01) %>%
    arrange(desc(logFC)) %>% .$symbol %>% as.character() %>% .[1:5]


#Genes to highlight
interest <- filter(topTab, symbol %in% c(topUp,topDown)) 
imminterest <- filter(topTab, symbol %in% immuno)
lacinterest <- filter(topTab, symbol %in% lac)

#VolcanoPlot
volcano <- ggplot(topTab,aes(x=logFC,y=-log10(FDR))) +
    geom_point(size=2,color="grey50",pch=20) +
    geom_hline(yintercept=2,lty="dashed") +
    geom_vline(xintercept=1,lty="dashed") +
    geom_vline(xintercept=-1,lty="dashed") +
    geom_point(data=interest, aes(x=logFC, y=-log10(FDR)), size=3, color="black",pch=20) +
    geom_point(data=lacinterest, aes(x=logFC, y=-log10(FDR)), size=3, color="dodgerblue",pch=20) +
    geom_point(data=imminterest, aes(x=logFC, y=-log10(FDR)), size=3, color="coral",pch=20) +
    geom_label_repel(data=interest, aes(x=logFC,y=-log10(FDR),label=symbol)) +
    xlab("Log2 Fold Change") +
    ylab("-Log10(P value)") 

#dummy legend
dum <- data.frame(x=c(1,1),y=c(1,1),grp=c("Immune response","Lactation"))
dumleg <- ggplot(dum, aes(x,y,color=grp)) +
    geom_point() +
    scale_color_manual(values=c("coral","dodgerblue")) +
    theme(legend.position="bottom",
	  legend.title=element_blank(),
	  legend.direction="horizontal") +
    guides(colour = guide_legend(override.aes = list(size=3))) 
dumleg <- get_legend(dumleg)


# ---- DEmaintainedInPI ----

#Normalize
m.norm <- t(t(m)/pD$sf)
rownames(m.norm) <- as.character(fD$symbol)
genes <- c("Csn2","Csn1s1","Csn1s2a","Csn3")

#Create DF to distinguish between 4-L and 4-PI
forPl <- data.frame(t(m.norm)[,c(genes)]+1,
	barcode=colnames(m.norm))
add <- select(pD,barcode, cluster, Condition) %>%
    mutate(barcode=as.character(barcode))
forPl <- left_join(forPl,add,by="barcode") %>%
    filter(cluster %in% c(5,4))
forPl <- melt(forPl,id=c("barcode","cluster","Condition")) %>%
    mutate(group=paste0(Condition,cluster)) %>%
    mutate(group=gsub("NP5|G5|L5|PI5","5",group)) %>%
    mutate(group=gsub("PI4","4-PI",group)) %>%
    mutate(group=gsub("L4","4-L",group)) %>%
    mutate(group=factor(group,levels=c("5","4-L","4-PI")))

#Plot
pal <- brewer.pal(n=9,name="Paired")[c(4,5)]
ExpPlot <- ggplot(forPl, aes(y=value,x=group,color=cluster)) +
    geom_jitter(size=0.9) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
		 geom = "crossbar", width = 1,color="black") +
    facet_grid(~variable) +
    ylab("Expression") +
    xlab("") +
    theme(strip.background=element_blank(),
	  strip.text=element_text(face="bold"),
	  legend.position="bottom",
	  legend.direction="horizontal",
	  ) +
    scale_colour_manual(values=pal)+
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_y_log10(breaks=c(10,25,50,100,200,300,500,1000))


# Combine plots
leg <- plot_grid(dumleg,get_legend(ExpPlot),nrow=1)
ExpPlot <- ExpPlot %+% guides(color=FALSE)
subp <- plot_grid(volcano,ExpPlot,nrow=1,labels=c("c","d"))
subp <- plot_grid(subp,leg,rel_heights=c(1,0.1),nrow=2)

fullP <- plot_grid(subP0,subp,nrow=2)

# cairo_pdf("../paper/figures/Figure4.pdf",width=10.75,height=15.19)
plot_grid(fullP,NULL,nrow=2,rel_heights=c(1,0.75))
# dev.off()

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


comps <- c(5,4)
out <- list()
library(edgeR)
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
    xxx <- topTags(res,n=Inf,sort.by="PValue")
    topTab <- xxx$table
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

library(cowplot)

interest <- filter(forPlot, Gene %in% c("Kit","Hey1","Cd14",
					"Prlr","Esr1","Pgr"))
p <- ggplot(forPlot, aes(x=NullParFC,y=ParousFC)) +
    geom_point(color="grey50",size=2) +
    #     geom_point(data=interest, aes(x=NullParFC, y=ParousFC), size=4, color="red",pch=1) +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), color="black") +
    geom_label_repel(data=filter(interest,NullParFC < 0), aes(x=NullParFC,y=ParousFC,label=Gene)) +
    geom_label_repel(data=filter(interest,NullParFC > 0), aes(x=NullParFC,y=ParousFC,label=Gene)) +
    xlab("LFC of C5 vs. luminal cells") +
    ylab("LFC of C4 vs. luminal cells") +
    geom_hline(yintercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    coord_equal(xlim=c(-5,5),ylim=c(-5,5))


subP0 <- plot_grid(simMat[[4]],p,labels=c("a","b"))


####
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

#########################################
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
    xxx <- topTags(res,n=Inf,sort.by="PValue")
    topTab <- xxx$table

lac <- c("Btn1a1","Lalba","B4galt1","Csn3","Csn1s2a","Csn1s1","Csn2","Hk2","Xdh","Vegfa")

immuno <- c("Hp","Slpi","H2-K1", "B2m", "H2-Q7", "Lbp", "Tlr2", "Ltf", "Ifit1",
	    "Cd1d1")

forxls <- select(topTab, id, symbol, logFC, unshrunk.logFC, logCPM, PValue, FDR)
write.csv(forxls,file="../paper/supps/DE_C4vsC5.csv",quote=FALSE)

topUp <- filter(topTab, FDR < 0.01) %>%
    arrange(logFC) %>% .$symbol %>% as.character() %>% .[1:5]
topDown <- filter(topTab, FDR < 0.01) %>%
    arrange(desc(logFC)) %>% .$symbol %>% as.character() %>% .[1:5]


interst <- filter(topTab, symbol %in% c(topUp,topDown)) 
imminterest <- filter(topTab, symbol %in% immuno)
lacinterest <- filter(topTab, symbol %in% lac)

volcano <- ggplot(topTab,aes(x=logFC,y=-log10(FDR))) +
    geom_point(size=2,color="grey50",pch=20) +
    geom_hline(yintercept=2,lty="dashed") +
    geom_vline(xintercept=1,lty="dashed") +
    geom_vline(xintercept=-1,lty="dashed") +
    geom_point(data=interst, aes(x=logFC, y=-log10(FDR)), size=3, color="black",pch=20) +
    geom_point(data=lacinterest, aes(x=logFC, y=-log10(FDR)), size=3, color="dodgerblue",pch=20) +
    geom_point(data=imminterest, aes(x=logFC, y=-log10(FDR)), size=3, color="coral",pch=20) +
    geom_label_repel(data=interst, aes(x=logFC,y=-log10(FDR),label=symbol)) +
    xlab("Log2 Fold Change") +
    ylab("-Log10(P value)") 

dum <- data.frame(x=c(1,1),y=c(1,1),grp=c("Immune response","Lactation"))
dumleg <- ggplot(dum, aes(x,y,color=grp)) +
    geom_point() +
    scale_color_manual(values=c("coral","dodgerblue")) +
    theme(legend.position="bottom",
	  legend.title=element_blank(),
	  legend.direction="horizontal") +
    guides(colour = guide_legend(override.aes = list(size=3))) 

dumleg <- get_legend(dumleg)


#####################################
#Normalize
m.norm <- t(t(m)/pD$sf)

stopifnot(identical(rownames(m.norm),as.character(fD$id)))
rownames(m.norm) <- as.character(fD$symbol)
###################################################
###################################################
genes <- c("Csn2","Csn1s1","Csn1s2a","Csn3")
forPl <- data.frame(t(m.norm)[,c(genes)]+1,
	barcode=colnames(m.norm))
add <- select(pD,barcode, cluster, Condition) %>%
    mutate(barcode=as.character(barcode))
forPl <- left_join(forPl,add,by="barcode") %>%
    filter(cluster %in% c(5,4))
forPl <- melt(forPl,id=c("barcode","cluster","Condition")) %>%
    mutate(group=paste0(Condition,cluster)) %>%
    mutate(group=gsub("V5|I5|L5|P5","5",group)) %>%
    mutate(group=gsub("I4","4-PI",group)) %>%
    mutate(group=gsub("L4","4-L",group)) %>%
    mutate(group=factor(group,levels=c("5","4-L","4-PI")))

library(RColorBrewer)
pal <- brewer.pal(n=9,name="Paired")[c(4,5)]
ExpPlot <- ggplot(forPl, aes(y=value,x=group,color=cluster)) +
    geom_jitter(size=0.9) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
		 geom = "crossbar", width = 1,color="black") +
    facet_grid(~variable) +
    ylab("Log-Expression") +
    xlab("") +
    theme(#axis.line.x=element_blank(),
	  #           axis.text.x=element_blank(),
	  #           axis.ticks.x=element_blank(),
	  strip.background=element_blank(),
	  strip.text=element_text(face="bold"),
	  legend.position="bottom",
	  legend.direction="horizontal",
	  ) +
    scale_colour_manual(values=pal)+
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_y_log10()



leg <- plot_grid(dumleg,get_legend(ExpPlot),nrow=1)
ExpPlot <- ExpPlot %+% guides(color=FALSE)
subp <- plot_grid(volcano,ExpPlot,nrow=1,labels=c("c","d"))
subp <- plot_grid(subp,leg,rel_heights=c(1,0.1),nrow=2)

fullP <- plot_grid(subP0,subp,nrow=2)

cairo_pdf("Figure4.pdf",width=10.75,height=15.19)
plot_grid(fullP,NULL,nrow=2,rel_heights=c(1,0.75))
dev.off()

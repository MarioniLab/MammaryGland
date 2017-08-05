# Figure 2

library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(scatterplot3d)
library(gridGraphics)
library(gridExtra)
source("functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
dms <- read.csv("../data/Robjects/secondRun_2500/dm_all.csv")
pD <- dataList[[2]]
pD <- right_join(pD,dms,by="barcode")


# Plot for luminal and basal cells
set.seed(300)
pD <- pD[sample(1:nrow(pD),nrow(pD)),]
cols <- pD$SuperColor

scatterplot3d(x=pD[,"DC1"],
	      y=-pD[,"DC2"],
	      z=pD[,"DC3"],
	      color=cols,
	      pch=20,
	      angle=60,
              cex.symbols=1.5,
	      scale.y=0.5,
	      mar=c(5,3,-0.1,3)+0.1,
	      xlab="Component 1",
	      ylab="-Component 2",
	      zlab="-Component 3",
	      box=FALSE)

g <- grab_grob()
g <- grid.arrange(g)
dev.off()

# Create the alternative view inlet
inlet <- ggplot(pD, aes(DC1,DC2,color=SuperCluster)) +
    geom_point(size=2,pch=20) +
    scale_color_manual(values=levels(cols))+
    xlab("Component 1") +
    ylab("-Component 2") 

library(ggExtra)
ggMarginal(inlet, type="histogram",margins="x",fill="white",bins=30)
# cairo_pdf("../paper/figures/F2inlet.pdf",width=17.54,height=17.54)
inlet
# dev.off()

# Create legend for luminal only and luminal/basal plot
pal <- levels(cols)
# forcLeg <- filter(pD, !cluster %in% c(4))
clustLeg <- ggplot(pD, aes(x=tSNE1,y=tSNE2,color=SuperCluster)) +
    geom_point() +
    scale_color_manual(values=pal) +
    theme(legend.direction="horizontal",
	  legend.title=element_blank(),
	  legend.position="bottom")+
    guides(colour=guide_legend(nrow=1,
			       override.aes=list(size=3)))
clustLeg <- get_legend(clustLeg)

# ---- LuminalOnly ----
dms <- read.csv("../data/Robjects/secondRun_2500/dm_luminal.csv")
pD <- dataList[[2]]
pD <- right_join(pD,dms,by="barcode")

cols <- levels(pD$SuperColors)

# Luminal compartment colored by clusters
p.clust <- ggplot(pD, aes(x=DC1,y=DC2, color=SuperCluster)) +
    geom_point(size=2, pch=20) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(values=cols)+
    guides(colour=FALSE) +
    xlab("Component 1") +
    ylab("Component 2") 

p.clust <- plot_grid(NULL,p.clust,NULL,ncol=1,rel_heights=c(0.3,1,0.30))
    

# ---- GeneExpressionTrends ----

# load data again
m <- dataList[[1]]
fD <- dataList[[3]]
m <- m[,pD$barcode]

m.norm <- t(t(m)/pD$sf)

# Genes to plot for trends
genes <- c("Aldh1a3","Csn2","Glycam1","Pgr","Esr1")
rownames(m.norm) <- fD$symbol
exps <- log2(m.norm[genes,]+1)
exps <- t(exps/apply(exps,1,max))

# setup df
fPlot <- data.frame(exps,
		    barcode=colnames(m.norm))
fPlot <- full_join(fPlot, pD[,c("barcode","DC1","DC2")], by="barcode")

# color
pal <- colorRampPalette(brewer.pal(n=7,name="YlOrRd"))(200)

#  forloop for all genes
pltlist <- list()
for (gene in genes) {
    fPl <- arrange_(fPlot, gene)
    p <- ggplot(fPl, aes_string(x="DC1",y="DC2", color=gene)) +
	geom_point(size=2, pch=20) +
	#         scale_color_gradient(high="#D73027",low="#4575B4") +
	scale_color_gradientn(colours=pal) +
	ggtitle(gene) +
	theme_void() +
	theme(legend.position="bottom",
	      legend.direction="horizontal",
	      legend.title=element_blank())
    pltlist[[gene]] <- p
}

leg <- get_legend(pltlist[[1]])
pls <- lapply(pltlist,function(x){
	      x <- x %+% guides(color=FALSE) 
	      return(x)})


# Combine all plots

subPa <-plot_grid(g,clustLeg,nrow=2,labels=c("a"),rel_heights=c(1,0.1))
subPb1 <- plot_grid(plotlist=pls,nrow=3)
subPb1b <- plot_grid(subPb1,leg,nrow=2,rel_heights=c(1,0.05))
subPb2 <- plot_grid(p.clust,subPb1b,labels=c("b"))

cairo_pdf("../paper/figures/Figure2.pdf",width=8.27,height=11.69)
plot_grid(subPa,subPb2,nrow=2)
dev.off()

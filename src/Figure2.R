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
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
dms <- read.csv("../data/Robjects/dm_all.csv")
pD <- dataList[[2]]
pD <- right_join(pD,dms,by="barcode")

# Set color scale according to F1a
pal <- brewer.pal(n=9,name="Paired")
cols <- mapvalues(pD$cluster,levels(pD$cluster)[-c(1,11)],
		  pal)

# Plot for luminal and basal cells
scatterplot3d(x=pD[,"DC1"],
	      y=-pD[,"DC2"],
	      z=-pD[,"DC3"],
	      color=cols,
	      pch=20,
	      angle=40,
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
pD$cluster <- factor(pD$cluster,levels=c(1,2,3,4,5,6,7,8,9))
pal2 <- pal[-4]
inlet <- ggplot(pD, aes(-DC3,DC1,color=cluster)) +
    geom_point(size=5) +
    xlab("Component 1") +
    ylab("-Component 3") +
    scale_color_manual(values=pal2) +
    guides(color=FALSE)

# cairo_pdf("../paper/figures/F2inlet.pdf",width=17.54,height=17.54)
inlet
# dev.off()

# Create legend for luminal only and luminal/basal plot
pal <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,6,7,8,9)]
forcLeg <- filter(pD, !cluster %in% c(4))
clustLeg <- ggplot(forcLeg, aes(x=tSNE1,y=tSNE2,color=cluster)) +
    geom_point() +
    scale_color_manual(values=pal) +
    theme(legend.direction="horizontal",
	  legend.title=element_blank(),
	  legend.position="bottom")+
    guides(colour=guide_legend(nrow=1,
			       override.aes=list(size=3)))
clustLeg <- get_legend(clustLeg)

# ---- LuminalOnly ----
dms <- read.csv("../data/Robjects/dm_luminal.csv")
pD <- dataList[[2]]
pD <- right_join(pD,dms,by="barcode")

clusts <- as.numeric(as.character(sort(unique(pD$cluster))))
cols <- brewer.pal(n=9,name="Paired")[clusts]
names(cols) <- clusts

# Luminal compartment colored by clusters
p.clust <- ggplot(pD, aes(x=DC1,y=DC2, color=cluster)) +
    geom_point(size=2, pch=20) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(values=cols)+
    guides(colour=FALSE) +
    xlab("Component 1") +
    ylab("Component 2") 

    

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

# cairo_pdf("../paper/figures/Figure2.pdf",width=8.27,height=11.69)
plot_grid(subPa,subPb2,nrow=2)
# dev.off()

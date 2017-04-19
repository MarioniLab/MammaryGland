# Figure 2

library(plyr)
library(splines)
library(destiny)
library(dplyr)
library(reshape)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(scatterplot3d)
library(gridGraphics)
library(gridExtra)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# Remove QC-fails,outlier and immune cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier 
m <- m[,keepCells]
pD <- pD[keepCells,]

# ---- BasalAndLuminalNPandG ----

condComb <- c("NP","G")
keepCells <- pD$Condition %in% condComb
m.vp <- m[,keepCells]
pD.vp <- pD[keepCells,]

# Remove lowly expressed genes
keep <- rowMeans(m.vp)>0.1
m.vp <- m.vp[keep,]
fD.vp <- fD[keep,]

# Normalize 
m.norm <- t(t(m.vp)/pD.vp$sf)

# Compute HVGs
brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
fD.vp$highVar <- fD.vp$id %in% brennecke

# Prepare expression matrix
exps <- m.norm[fD.vp$highVar,]
exps <- t(log(exps+1))

# Compute diffusion map
set.seed(rnd_seed)
dm <- DiffusionMap(exps,n_eigs=20,k=50)
dms <- eigenvectors(dm)[,1:3]

# Set color scale according to F1a
pal <- brewer.pal(n=9,name="Paired")
cols <- mapvalues(pD.vp$cluster,levels(pD.vp$cluster)[-c(1,10)],
		  pal)

# Plot for luminal and basal cells
scatterplot3d(x=dms[,1],
	      y=-dms[,2],
	      z=-dms[,3],
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

#Create the alternative view inlet
inletD <- pD.vp
inletD$dc1 <- dms[,1]
inletD$dc3 <- dms[,3]
inletD$cluster <- factor(inletD$cluster,levels=c(1,2,3,4,5,6,7,8,9))
pal2 <- pal[-4]
inlet <- ggplot(inletD, aes(-dc3,dc1,color=cluster)) +
    geom_point(size=5) +
    xlab("Component 1") +
    ylab("-Component 3") +
    scale_color_manual(values=pal2) +
    guides(color=FALSE)

cairo_pdf("../paper/figures/F2inlet.pdf",width=17.54,height=17.54)
inlet
dev.off()

# ---- LuminalOnlyNPandG ----

# Cells
keepCells <- !(pD$cluster %in% c(6,7,9)) & pD$Condition %in% condComb
m.vp <- m[,keepCells]
pD.vp <- pD[keepCells,]

# Genes 
keep <- rowMeans(m.vp)>0.1
m.vp <- m.vp[keep,]
fD.vp <- fD[keep,]

# Normalize
m.norm <- t(t(m.vp)/pD.vp$sf)

# HVGs
brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
fD.vp$highVar <- fD.vp$id %in% brennecke

# Prepare expression matrix by selecting only HVGs and log-transformation
exps <- m.norm[fD.vp$highVar,]
exps <- t(log(exps+1))

# Compute Diffusion map
set.seed(rnd_seed)
dm <- DiffusionMap(exps,n_eigs=20,k=50)
dcs <- eigenvectors(dm)[,1:2]
pD.vp$DC1 <- eigenvectors(dm)[,1]
pD.vp$DC2 <- eigenvectors(dm)[,2]

# Color scheme for plots
clusts <- as.numeric(as.character(sort(unique(pD.vp$cluster))))
cols <- brewer.pal(n=9,name="Paired")[clusts]
names(cols) <- clusts

# Luminal compartment colored by clusters
p.clust <- ggplot(pD.vp, aes(x=DC1,y=DC2, color=cluster)) +
    geom_point(size=2, pch=20) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(values=cols)+
    guides(colour=FALSE) +
    xlab("Component 1") +
    ylab("Component 2") 

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
    

# ---- GeneExpressionTrends ----

genes <- c("Aldh1a3","Kit","Csn2","Glycam1","Pgr","Esr1")
rownames(m.norm) <- fD.vp$symbol
exps <- log2(m.norm[genes,]+1)
exps <- t(exps/apply(exps,1,max))
fPlot <- data.frame(exps,
		    barcode=colnames(m.norm))
fPlot <- join(fPlot, pD.vp[,c("barcode","DC1","DC2")], by="barcode")
pal <- colorRampPalette(brewer.pal(n=7,name="YlOrRd"))(200)
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

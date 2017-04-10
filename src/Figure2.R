#########################################
#
#Figure 2
#
#########################################

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
p <- ggplot(pD.vp, aes(x=DC1,y=DC2, color=cluster)) +
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

# Aldehyde trend along component2
ald <- filter(fD.vp, symbol %in% c("Aldh1a3")) %>% .$id
fDC2Trend <- data.frame("DC2"=pD.vp$DC2,
		   "Aldh1a3"=unname(log2(t(m.norm)[,ald]+1)))
DC2Trend <- ggplot(fDC2Trend, aes(x=DC2, y=Aldh1a3)) +
    geom_point(size=0.5,color="#7570b3",alpha=0.5) +
    geom_smooth(method="glm",formula="y~ns(x,df=5)",
		method.args=list(family="quasipoisson"),color="#7570b3") +
    ylab("Log-Expression") +
    xlab("Component 2") +
    coord_flip()

# Csn2 and Pgr trend along component1
genes <- filter(fD.vp, symbol %in% c("Csn2","Pgr")) %>% .$id
fDC1Trend <- data.frame("DC1"=pD.vp$DC1,
		   "Csn2"=unname(log2(t(m.norm)[,genes[1]]+1)),
		   "Pgr"=unname(log2(t(m.norm)[,genes[2]]+1)),
		   "Aldh1a3"=NA) # include Aldh1a3 as NA to have a common legend
fDC1Trend <- melt(fDC1Trend,id="DC1",variable_name="Gene")
DC1Trend <- ggplot(fDC1Trend, aes(x=DC1, y=value, color=Gene)) +
    geom_point(size=0.5,alpha=0.5) +
    geom_smooth(method="glm",formula="y~ns(x,df=5)",
		method.args=list(family="quasipoisson")) +
    scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")) +
    ylab("Log-Expression") +
    xlab("Component 1") +
    guides(color=guide_legend(override.aes=list(fill=NA,size=2))) +
    theme(legend.direction="vertical",legend.title=element_blank()) +
    xlab("") 

# Combine all plots
leg2 <- get_legend(DC1Trend)

DC1Trend <- DC1Trend %+% guides(colour=FALSE)

subP <- plot_grid(p,DC2Trend,DC1Trend,align="hv",rel_widths=c(1,0.6),
		  rel_heights=c(1,0.6))

subP1 <- subP + draw_grob(leg2,0.7,-0.11,1/3,0.5)
subP1 <- plot_grid(subP1,labels=c("b"))

subPab <-plot_grid(g,labels=c("a"))
subP0 <- plot_grid(subPab,clustLeg,nrow=2,rel_heights=c(1,0.1))

cairo_pdf("../paper/figures/Figure2.pdf",width=8.41,height=12.54)
plot_grid(subP0,subP1,labels=c("",""),nrow=2)
dev.off()

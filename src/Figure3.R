# Figure 3

library(plyr)
library(cowplot)
library(pheatmap)
library(splines)
library(dplyr)
library(destiny)
library(ggplot2)
source("functions.R")
library(RColorBrewer)
library(viridis)
library(lmtest)
library(gridExtra)

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

dms <- read.csv("../data/Robjects/dm_luminal.csv")
pD <- right_join(pD,dms,by="barcode")

# ---- PlotWithBranchDPT -----

b1 <- pD[pD$branch %in% c("Root","Intermediate",
			  "Hormone-sensing lineage"),]
pb1 <- ggplot(pD,aes(x=DC1,y=-DC2)) +
    geom_point(size=0.8,color="grey80") +
    geom_point(data=b1,aes(x=DC1,y=-DC2,color=dpt)) +
    scale_color_viridis(option="magma",begin=1,end=0) +
    xlab("Component 1") +
    ylab("-Component 2") +
    ggtitle("Hormone-sensing lineage") +
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank(),
	  legend.title=element_blank())


b2 <- pD[pD$branch %in% c("Root","Intermediate",
			  "Secretory lineage"),]
pb2 <- ggplot(pD,aes(x=DC1,y=-DC2)) +
    geom_point(size=0.8,color="grey80") +
    geom_point(data=b2,aes(x=DC1,y=-DC2,color=dpt)) +
    scale_color_viridis(option="magma",begin=1,end=0) +
    xlab("Component 1") +
    ylab("-Component 2") +
    ggtitle("Secretory lineage") +
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank(),
	  legend.title=element_blank())

branches <- plot_grid(pb2,NULL,pb1,NULL,nrow=1,
		      rel_widths=c(1,.2,1,.4))

# ---- BranchDE ----

out <- readRDS("../data/Robjects/BranchDEList.rds")

# ---- Heatmap -----

m.hrm <- out[[1]][["mSmooth"]]
m.alv <- out[[2]][["mSmooth"]]

# Combine smoothed expression values for heatmap
m.both <- cbind(m.alv[,c(ncol(m.alv):1)],m.hrm) # reverse order of alveolar cells for heatmap

# Scale expression values
m.both <- t(scale(t(m.both)))
m.both[m.both>3] <- 3 # cut at 3 for visualization
m.both[m.both<-3] <- -3 # cut at -3 for visualization


# Matrix for heatmap
genes <- c(out[["genes.sameGrad"]][1:50],
	   out[["genes.diffGrad"]][1:50])

m.heat <- m.both[genes,]

# Annotation Data
pD.ord <- pD
rownames(pD.ord) <- pD.ord$barcode
pD.ord <- pD.ord[colnames(m.heat),]
annoCol <- data.frame("Cluster"=pD.ord$cluster,
		      "Pseudotime"=rank(pD.ord$dpt,ties.method="first")
		      )

# Rename the alveolar cells, so that there are no duplicate names for rows
colnames(m.heat) <- c(paste0("Alv.",colnames(m.heat)[1:ncol(m.alv)]),
		       colnames(m.heat)[(ncol(m.alv)+1):ncol(m.heat)])
rownames(annoCol) <- colnames(m.heat)

#Set colorscheme for heatmap
clustCol <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,8)]
names(clustCol) <- c(1,2,3,5,8)
dptcols <- viridis(n=nrow(annoCol),,option="magma",begin=1,end=0)
names(dptcols) <- c(1:length(dptcols))
annoColors <- list("Cluster"=clustCol,
		   "Pseudotime"=dptcols)

#Plot heatmap
p0 <- pheatmap(m.heat,
	 cluster_cols=FALSE,
	 cluster_rows=FALSE,
	 annotation=annoCol,
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 annotation_legend=FALSE,
	 gaps_col=ncol(m.alv),
	 gaps_row=c(50),
	 show_rownames=TRUE,
	 fontsize=6)


# ---- PlotTFExamples ----


#Set Genes to plot (TFs that are DE)
features <- c("Creb5","Hey1","Fosl1",
	      "Runx1","Tox2","Bhlhe41",
	      "Elf5","Foxs1","Ehf")

#extract expression for first branch
p1 <- out[[1]][["pD"]] %>% arrange(DPTRank)
yhet1 <- data.frame(t(m.hrm)[,features])
yhet1$barcode <- as.character(p1$barcode)
raw1 <- data.frame(t(out[[1]][["m"]][features,yhet1$barcode]))
colnames(raw1) <- paste0("raw",colnames(raw1))
raw1$barcode <- yhet1$barcode
fplot1 <- join(p1,yhet1,by="barcode") 
fplot1 <- join(fplot1,raw1,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))

#extract expression for second branch
p2 <- out[[2]][["pD"]] %>% arrange(DPTRank)
yhet2 <- data.frame(t(m.alv)[,features])
yhet2$barcode <- as.character(p2$barcode)
raw2 <- data.frame(t(out[[2]][["m"]][features,yhet2$barcode]))
colnames(raw2) <- paste0("raw",colnames(raw2))
raw2$barcode <- yhet2$barcode
fplot2 <- join(p2,yhet2,by="barcode") 
fplot2 <- join(fplot2,raw2,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))


#Plot dpt-dependent expression for features
pList <- list()
for (feature in features) {
    pnts <- paste0("raw",feature)
    lns <- feature
    p <- ggplot() +
	geom_point(size=0.8,data=fplot1,aes_string(x="dptNorm",y=pnts, color="cluster")) +
	geom_point(size=0.8,data=fplot2,aes_string(x="dptNorm",y=pnts, color="cluster")) +
	geom_line(data=fplot1,aes_string(x="dptNorm",y=lns),lty="dashed") +
	geom_line(data=fplot2,aes_string(x="dptNorm",y=lns)) +
	ggtitle(feature) +
	ylab("") +
	scale_color_manual(values=clustCol) +
	scale_x_continuous(breaks=c(0,1)) +
	theme(legend.position="bottom",
	      legend.direction="horizontal",
	      legend.title=element_blank()) +
	guides(colour = guide_legend(override.aes = list(size=3))) +
	xlab("") 
    pList[[feature]] <- p
}

# create a dummy legend for dasehd/solid line
dummyD <- data.frame(x=rnorm(10),y=rnorm(10),
		     branch=c(rep("Secretory lineage",5),rep("Hormone-sensing lineage",5)))
dummyP2 <- ggplot(dummyD,aes(x,y,lty=branch)) +
    geom_line() +
    theme(legend.position="bottom",
	  #           legend.direction="horizontal",
	  legend.title=element_blank()) +
    scale_linetype_manual(values=c("dashed","solid"))
legs1 <- get_legend(dummyP2)


# extract color legend from plotList and delete legend
legs <- get_legend(pList[[1]])
legs <- plot_grid(legs1,legs,nrow=1)
pls <- lapply(pList,function(x){
	      x <- x %+% guides(color=FALSE) 
	      return(x)})

# add Xlab/Ylab to bottom/left center plots 
pls[["Foxs1"]] <- pls[["Foxs1"]] %+% xlab("Normalized Pseudotime")
pls[["Runx1"]] <- pls[["Runx1"]] %+% ylab("Log-Expression")

#Combine all plots in a single figure
expPlot <- plot_grid(plotlist=pls,ncol=3,labels=c("c","","",
						  "d","","",
						  "e","",""))

expPlot <- plot_grid(expPlot,legs,rel_heights=c(1,0.05),ncol=1)
expPlot <- plot_grid(NULL,expPlot,ncol=1,rel_heights=c(0.3,1))
htmp <- plot_grid(NULL,p0[[4]],ncol=1,rel_heights=c(.275,1),labels=c("a","b"),
		  vjust=c(1.5,0.5))

#Draw branches above heatmap
pb1 <- grid.arrange(pb1)
pb2 <- grid.arrange(pb2)
htmp <- htmp + draw_grob(pb2,-0.02,0.78,.35,.2)
htmp <- htmp + draw_grob(pb1,0.4,0.78,.35,.2)

fullP <- plot_grid(htmp,expPlot,ncol=2,rel_widths=c(1,1))

#close graphics device before plotting
# dev.off()
# cairo_pdf("../paper/figures/Figure3.pdf",width=16.55,height=13.0575)
fullP
# dev.off()
# 

# S4
library(pheatmap)
library(dplyr)
library(viridis)
library(cowplot)
library(RColorBrewer)

out <- readRDS("../data/Robjects/secondRun_2500/BranchDEList.rds")
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_Clustered.rds")
dms <- read.csv("../data/Robjects/secondRun_2500/dm_luminal.csv")

pD <- dataList[[2]]
pD <- right_join(pD,dms,by="barcode")


m.hrm <- out[[1]][["mSmooth"]]
m.alv <- out[[2]][["mSmooth"]]

# Combine smoothed expression values for heatmap
m.both <- cbind(m.alv[,c(ncol(m.alv):1)],m.hrm) # reverse order of alveolar cells for heatmap

# Scale expression values
m.both <- t(scale(t(m.both)))
m.both[m.both>3] <- 3 # cut at 3 for visualization
m.both[m.both<-3] <- -3 # cut at -3 for visualization

# Annotation Data
pD.ord <- pD
rownames(pD.ord) <- pD.ord$barcode
pD.ord <- pD.ord[colnames(m.both),]
annoCol <- data.frame("Cluster"=pD.ord$cluster,
		      "Pseudotime"=rank(pD.ord$dpt,ties.method="first")
		      )
#Rename the alveolar cells, so that there are no duplicate names for rows
colnames(m.both) <- c(paste0("Alv.",colnames(m.both)[1:ncol(m.alv)]),
		       colnames(m.both)[(ncol(m.alv)+1):ncol(m.both)])
rownames(annoCol) <- colnames(m.both)

#Set colorscheme for heatmap
clustCol <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,8)]
names(clustCol) <- c(1,2,3,5,8)
dptcols <- viridis(n=nrow(annoCol),,option="magma",begin=1,end=0)
names(dptcols) <- c(1:length(dptcols))
annoColors <- list("Cluster"=clustCol,
		   "Pseudotime"=dptcols)


m.heat1 <- m.both[out[["genes.sameGrad"]],]
m.heat2 <- m.both[out[["genes.diffGrad"]],]

##Plot
p0 <- pheatmap(m.heat1,
	 cluster_cols=FALSE,
	 cluster_rows=TRUE,
	 clustering_distance_rows="euclidean",
	 annotation=annoCol,
	 clustering_method="ward.D2",
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 treeheight_row=0,
	 legend=FALSE,
	 annotation_legend=FALSE,
	 gaps_col=ncol(m.alv),
	 show_rownames=FALSE)

p1 <- pheatmap(m.heat2,
	 cluster_cols=FALSE,
	 cluster_rows=TRUE,
	 clustering_distance_rows="euclidean",
	 annotation=annoCol,
	 clustering_method="ward.D2",
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 treeheight_row=0,
	 annotation_legend=TRUE,
	 gaps_col=ncol(m.alv),
	 show_rownames=FALSE)

# png("../paper/figures/S4.png",height=1248,width=1248)
comb <- plot_grid(p0[[4]],NULL,p1[[4]],nrow=1,vjust=0.5,rel_widths=c(1,.1,1))
comb
# dev.off()

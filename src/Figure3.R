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

dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

dms <- read.csv("../data/Robjects/secondRun_2500/dm_luminal.csv")
pD <- right_join(pD,dms,by="barcode")

# ---- PlotWithBranchDPT -----

b1 <- pD[pD$branch %in% c("Root","Intermediate",
			  "Hormone-sensing lineage"),]
pb1 <- ggplot(pD,aes(x=DC1,y=DC2)) +
    geom_point(size=0.8,color="grey80") +
    geom_point(data=b1,aes(x=DC1,y=DC2,color=dpt)) +
    scale_color_viridis(option="magma",begin=1,end=0) +
    xlab("Component 1") +
    ylab("-Component 2") +
    ggtitle("Hormone-sensing lineage") +
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank(),
	  legend.title=element_blank())


b2 <- pD[pD$branch %in% c("Root","Intermediate",
			  "Secretory lineage"),]
pb2 <- ggplot(pD,aes(x=DC1,y=DC2)) +
    geom_point(size=0.8,color="grey80") +
    geom_point(data=b2,aes(x=DC1,y=DC2,color=dpt)) +
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

out <- readRDS("../data/Robjects/secondRun_2500/BranchDEList.rds")

# ---- Heatmap -----
res <- list()
for (i in c(1,2)) {
    # Scale expression values
    m.heat <- t(scale(t(out[[i]][["mSmooth"]])))
    m.heat[m.heat>3] <- 3 # cut at 3 for visualization
    m.heat[m.heat<-3] <- -3 # cut at -3 for visualization


    # Matrix for heatmap
    x <- out[[i]][["Results"]]
    genes <- intersect(as.character(x[x$PAdjust < 0.01,"Gene"]),as.character(out[["genes.diffGrad"]]))
    m.heat <- m.heat[genes,]

    # Annotation Data
    pD.ord <- pD
    rownames(pD.ord) <- pD.ord$barcode
    pD.ord <- pD.ord[colnames(m.heat),]
    annoCol <- data.frame("Cluster"=factor(pD.ord$SubCluster),
			  "Pseudotime"=rank(pD.ord$dpt,ties.method="first")
			  )
    rownames(annoCol) <- rownames(pD.ord)

    #Set colorscheme for heatmap
    clustCol <- levels(factor(pD.ord$Colors))
    names(clustCol) <- levels(factor(pD.ord$SubCluster))
    dptcols <- viridis(n=nrow(annoCol),,option="magma",begin=1,end=0)
    names(dptcols) <- c(1:length(dptcols))
    annoColors <- list("Cluster"=clustCol,
		       "Pseudotime"=dptcols)
    if(i==1) {
	k <- 3
    } else {
	k <- 4
    }

    # clustering
    dis <- dist(m.heat)
    library(dynamicTreeCut)
    clus <- dynamicCluster(dis, lk="ward.D2", output="cluster")
    #     m.heat <- m.heat[clus$tree$order,]
    clusOrder <- unique(unname(clus$cluster[clus$tree$order]))
    rowGaps <- cumsum(unname(table(clus$cluster)[clusOrder]))

    
    m.heat <- m.heat[clus$tree$order,]
    #Plot heatmap
    p0 <- pheatmap(m.heat,
	     cluster_cols=FALSE,
	     cluster_rows=FALSE,
	     #              clustering_distance_rows="euclidean",
	     #              clustering_method="ward.D2",
	     annotation=annoCol,
	     show_colnames=FALSE,
	     #              cutree_rows=k,
	     gaps_row=rowGaps,
	     annotation_colors=annoColors,
	     annotation_legend=FALSE,
	     show_rownames=FALSE,
	     fontsize=6)

    gene.clusters <- clus$cluster
    names(gene.clusters) <- clus$tree$labels
    res[[i]] <- list(p0,gene.clusters)
}

# ---- Save list for GO-Analysis ----
hrm.genes <- data.frame("Branch"="Hrm",
		  "Cluster"=res[[1]][[2]],
		  "Gene"=names(res[[1]][[2]]))

alv.genes <- data.frame("Branch"="Alv",
		  "Cluster"=res[[2]][[2]],
		  "Gene"=names(res[[2]][[2]]))

write.csv(rbind(hrm.genes,alv.genes),file="../data/Robjects/secondRun_2500/BranchDECluster.R",
	  row.names=FALSE)


# ---- PlotTFExamples ----

#Set Genes to plot (TFs that are DE)
features <- c("Tead1","Fosl1","Hmga1",
	      "Runx1","Tox2","Bhlhe41",
	      "Elf5","Foxs1","Ehf")
# both matrices
m.hrm <- out[[1]][["mSmooth"]]
m.alv <- out[[2]][["mSmooth"]]

#extract expression for first branch
p1 <- out[[1]][["pD"]] %>% arrange(DPTRank)
yhet1 <- data.frame(t(m.hrm)[,features],
		    check.names=FALSE)
yhet1$barcode <- as.character(p1$barcode)
raw1 <- data.frame(t(out[[1]][["m"]][features,yhet1$barcode]),
		   check.names=FALSE)
colnames(raw1) <- paste0("raw",colnames(raw1))
raw1$barcode <- yhet1$barcode
fplot1 <- join(p1,yhet1,by="barcode") 
fplot1 <- join(fplot1,raw1,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))

#extract expression for second branch
p2 <- out[[2]][["pD"]] %>% arrange(DPTRank)
yhet2 <- data.frame(t(m.alv)[,features],
		    check.names=FALSE)
yhet2$barcode <- as.character(p2$barcode)
raw2 <- data.frame(t(out[[2]][["m"]][features,yhet2$barcode]),
		   check.names=FALSE)
colnames(raw2) <- paste0("raw",colnames(raw2))
raw2$barcode <- yhet2$barcode
fplot2 <- join(p2,yhet2,by="barcode") 
fplot2 <- join(fplot2,raw2,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))


#set colorscale
clustCol <- levels(fplot2$Colors)[levels(fplot2$SubCluster)[c(-12,-17,-18,-19,-20)] %in% unique(c(as.character(fplot2$SubCluster),as.character(fplot1$SubCluster)))]

fplot3 <- rbind(fplot1,fplot2)
#Plot dpt-dependent expression for features
pList <- list()
for (feature in features) {
    pnts <- ggname(paste0("raw",feature))
    lns <- ggname(feature)
    p <- ggplot() +
	geom_point(size=0.8,data=fplot3,aes_string(x="dptNorm",y=pnts, color="SubCluster")) +
	geom_line(data=fplot1,aes_string(x="dptNorm",y=lns),lty="dashed") +
	geom_line(data=fplot2,aes_string(x="dptNorm",y=lns)) +
	ggtitle(feature) +
	ylab("") +
	scale_color_manual(values=clustCol) +
	scale_x_continuous(breaks=c(0,1)) +
	theme(legend.position="bottom",
	      legend.direction="horizontal",
	      legend.title=element_blank(),
	      plot.title=element_text(face="italic")) +
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
expPlot <- plot_grid(plotlist=pls,ncol=3,labels=c("b","","",
						  "c","","",
						  "d","",""))
expPlot <- plot_grid(expPlot,legs,rel_heights=c(1,0.05),ncol=1)
branches <- plot_grid(NULL,pb1,pb2,NULL,rel_widths=c(0.3,1,1,0.3),nrow=1,labels=c("a"))
expPlot <- plot_grid(branches,expPlot,ncol=1,rel_heights=c(0.3,1))

htmps <- plot_grid(res[[1]][[1]][[4]],
		   NULL,
		   res[[2]][[1]][[4]],
		   NULL,
		   nrow=2,
		   scale=0.9,
		   labels=c("e","","f",""),
		   rel_widths=c(1,0.75,1,0.75))
fullP <- plot_grid(expPlot,htmps,ncol=2)


#close graphics device before plotting
dev.off()
cairo_pdf("../paper/figures/Figure3.pdf",width=16.55,height=13.0575)
fullP
dev.off()
# 

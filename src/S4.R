library(scran)
library(plyr) 
library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
library(cowplot)
library(dynamicTreeCut)
library(Rtsne)
library(pheatmap)
source("functions.R")

rnd_seed <- 300
clustDat <- readRDS("../data/Robjects/ClusterComparison.rds") 
bootDat <- readRDS("../data/Robjects/ClusterBootstrap.rds")
dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")

#Commence
m.full <- dataList[[1]]
pD.full <- dataList[[2]]
fD.full <- dataList[[3]]
#Gene and Cell filtering
m <- m.full[fD.full$keep,pD.full$PassAll]
pD <- pD.full[pD.full$PassAll,]
fD <- fD.full[fD.full$keep,]
#relevel pD Condition to get it in logical sequence
pD$Condition <- mapvalues(pD$Condition, from=c("V","P","L","I"),
			  to=c("NP", "G",
			       "L", "PI"))
pD$Condition <- factor(pD$Condition, levels=c("NP","G","L","PI"))
pD$SampleID <- mapvalues(pD$SampleID, from=c("V1","V2","P1","P2","L1","L2","I1","I2"),
			  to=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2"))
pD$SampleID <- factor(pD$SampleID,levels=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2")) 

#Normalization
m <- t(t(m)/pD$sf)
m.sub <- m[fD$highVar,]

dm <- "euclidean"
lk <- "average"
ds <- 1
minSize <- 15
clust <- dynamicCluster(t(m.sub), dm=dm,
				 lk=lk,
				 ds=ds,
				 output="Cluster",
				 minSize=minSize)

## Extract as grob for plot
library(gridGraphics)
library(gridExtra)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

ds <- 0
# Bootstrap stability of clusters
slct <- paste(dm,lk,ds,sep="_")
boot <- bootDat[[slct]]$bootresult
boot <- boot[-1,]
bootx <- t(boot)
colnames(bootx) <- c("a","b","c","d","e","f","g","h")
boxplot(bootx,ylim=c(0,1))
title("Deep Split 0")
g1 <- grab_grob()

ds <- 1
# Bootstrap stability of clusters
slct <- paste(dm,lk,ds,sep="_")
boot <- bootDat[[slct]]$bootresult
boot <- boot[c(2,3,4,5,6,7,9,10,11,8),]
boxplot(t(boot),ylim=c(0,1))
title("Deep Split 1")
g2 <- grab_grob()

ds <- 2
# Bootstrap stability of clusters
slct <- paste(dm,lk,ds,sep="_")
boot <- bootDat[[slct]]$bootresult
bootx <- t(boot)
colnames(bootx) <- paste0("X",1:ncol(bootx))
boxplot(bootx,ylim=c(0,1))
title("Deep Split 2")
g3 <- grab_grob()


g <- grid.arrange(g1,g2,g3,ncol=3)

#Similarity matrix
pD$cluster <- clust$cluster
pD$cluster <- mapvalues(pD$cluster,c(0,1,2,3,4,5,6,7,9,10,8),
			c(0,1,2,3,4,5,6,7,8,9,10))
pD$cluster <- factor(pD$cluster,levels=c(0,1,2,3,4,5,6,7,8,9,10))
ord <- arrange(pD, cluster) 
ord$cluster <- as.numeric(as.character(ord$cluster))
simMat <- asSim(clust$dis)
simMat.ord <- as.matrix(simMat)[as.character(ord$barcode),
				as.character(ord$barcode)]
library(pheatmap)
annoCol <- data.frame("cluster"=as.factor(ord$cluster),
		      "Condition"=ord$Condition)
rownames(annoCol) <- as.character(ord$barcode)
library(RColorBrewer)
pal <- brewer.pal(n=11,"Paired")
pal <- pal[c(11,1,2,3,4,5,6,7,8,9,10)]
ann_colors <- list("cluster"=pal,
		   "Condition"=c("#F8766D","#7CAE00","#00BFC4","#C77Cff"))
names(ann_colors$cluster) <- c(min(ord$cluster):max(ord$cluster))
names(ann_colors$Condition) <- c("NP","G","L","PI")
# pheat <- pheatmap(simMat.ord,
#          cluster_rows=FALSE,
#          cluster_cols=FALSE,
#          show_rownames=FALSE,
#          show_colnames=FALSE,
#          annotation_col=annoCol,
#          annotation_colors=ann_colors)
# 
nclust <- length(unique(clust$cluster))
palette(pal)
htreeWithBars(clust,pD)
g4 <- grab_grob()
g4 <- grid.arrange(g4)
# 
# Plot of a few markers

genes <- c("Cd52","Cd74","Cd72")
m2 <- t(t(m.full[,pD.full$PassAll])/pD$sf)
p3 <- plotGeneDist(m2, pD, fD.full, genes, colorBy="cluster")
p3 <- p3 %+% facet_grid(variable~.) %+% xlab("Cluster") %+% ylab("Log-Expression")

out <- data.frame("C1"=numeric(nrow(m)))
for (clust in levels(pD$cluster)[-1]) {
    expr <- rowMeans(log2(m[,pD$cluster==clust]+1))
    colname <- paste0("C",clust)
    out[,colname] <- expr
}
meandist <- as.dist((1-cor(out))/2)
meandist <- dist(t(out))
x <- prcomp(t(out))
pal <- brewer.pal(n=10,"Paired")
palette(pal)
plot(x$x[,1],x$x[,2],col=c(1:10),pch=20,cex=2,xlab="PC1",ylab="PC2")
legend(x = 20, y = -10, legend = c(1:10), col = c(1:10), pch = 16)
g5 <- grab_grob()
g5 <- grid.arrange(g5)

pD <- filter(pD, cluster!=0)
LibrarySize <- ggplot(pD, aes(x=cluster,y=UmiSums)) +
    geom_violin(draw_quantiles=0.5)+
    ylab("Total Number of Molecules") +
    scale_y_log10() 

sub0 <- plot_grid(g,g4,ncol=1,labels=c("a","b"))
sub2 <- plot_grid(LibrarySize,p3,g5,nrow=1,labels=c("c","d","e"))

fin <- plot_grid(sub0,sub2,ncol=1,rel_heights=c(2,1))

cairo_pdf("../paper/figures/S3.pdf",height=14.28,width=14.28)
fin
dev.off()

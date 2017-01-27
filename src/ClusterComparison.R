library(scran)
library(dplyr)
library(knitr)
library(ggplot2)
library(dynamicTreeCut)
library(Rtsne)
library(pheatmap)
source("functions.R")

rnd_seed <- 300
clustDat <- readRDS("../data/Robjects/ClusterComparison.rds") 
bootDat <- readRDS("../data/Robjects/ClusterBootstrap.rds")
dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")


#Commence
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
#Gene and Cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]
#Normalization
clusters <- quickCluster(m)
pD$sf <- computeSumFactors(m,clusters=clusters)
m <- t(t(m)/pD$sf)
#Identification of HVG with various methods 
fD <- mutate(fD, meanExpression=rowMeans(m),
	     vars=apply(m,1,var),
	     CV2=apply(m,1,cv2),
	     dm=DM(meanExpression,CV2))
fit <- trendVar(log2(m+1))
decVar <- decomposeVar(log2(m+1),fit)
feat <- decVar[decVar$FDR < 0.1,] %>% na.omit() %>% rownames()
corfeat <- correlatePairs(m[feat,])
sig.feat <- corfeat$FDR <= 0.05
chosenfeat <- unique(c(corfeat$gene1[sig.feat],corfeat$gene2[sig.feat]))

# library(RBGL)
# g <- ftM2graphNEL(cbind(corfeat$gene1,corfeat$gene2)[sig.feat,], W=NULL, V=NULL, edgemode="undirected")
# cl <- highlyConnSG(g)$clusters
# cl <- cl[order(lengths(cl), decreasing=TRUE)]

highVar <- dplyr::arrange(fD,desc(dm)) %>% .$id
topX <- floor(0.05*length(highVar))

brennecke <- BrenneckeHVG(m,fdr=0.1,minBiolDisp=0.25)

fD$highDm <- fD$id %in% highVar[1:topX]
fD$highVar <- fD$id %in% chosenfeat
fD$brennecke <- fD$id %in% brennecke

genes2lbl <- c("Acta2", "Krt18", "Csn2", "Wap","Krt14")
lblDat <- filter(fD, symbol %in% genes2lbl)
V1 <- ggplot(arrange(fD,highDm), aes(x=meanExpression,y=CV2, color=highDm)) +
    geom_label(data=lblDat, aes(x=meanExpression, y=CV2, label=symbol), nudge_y=0.1) +
    geom_point(data=lblDat, aes(x=meanExpression, y=CV2), pch=1, size=4) +
    geom_point(alpha=0.8) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
V2 <- ggplot(arrange(fD,highVar), aes(x=meanExpression,y=CV2, color=highVar)) +
    geom_label(data=lblDat, aes(x=meanExpression, y=CV2, label=symbol), nudge_y=0.1) +
    geom_point(data=lblDat, aes(x=meanExpression, y=CV2), pch=1, size=4) +
    geom_point(alpha=0.8) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
V3 <- ggplot(arrange(fD,brennecke), aes(x=meanExpression,y=CV2, color=brennecke)) +
    geom_label(data=lblDat, aes(x=meanExpression, y=CV2, label=symbol), nudge_y=0.1) +
    geom_point(data=lblDat, aes(x=meanExpression, y=CV2), pch=1, size=4) +
    geom_point(alpha=0.8) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
library(cowplot)
plot_grid(V1,V2,V3,align="hv",ncol=1)

#Which Linkage is best
m.sub <- m[fD$brennecke,]
trafM <- log2(m.sub+1)
dms <- c("pearson","spearman","euclidean")
res <- NULL
for (dm in dms) {
    #Dissimilarity Measure
    if(dm!="euclidean") {
    dis <- as.dist((1-cor(trafM,method=dm))/2)
    } else {
	dis <- dist(t(trafM),method=dm)
    }

    #Hierarchical Clustering
    htree.average <- hclust(dis, method="average")
    htree.single <- hclust(dis, method="single")
    htree.complete <- hclust(dis, method="complete")
    htree.ward <- hclust(dis, method="ward.D2")

    ccc.average <- cor(cophenetic(htree.average),dis)
    ccc.single <- cor(cophenetic(htree.single),dis)
    ccc.complete <- cor(cophenetic(htree.complete),dis)
    ccc.ward <- cor(cophenetic(htree.ward),dis)

    tmp <- data.frame("Linkage"=c("Average",
					 "Single",
					 "Complete",
					 "Ward.D2"),
			     "Cophenetic Correlation Coefficient"=
				 c(ccc.average,ccc.single,
				   ccc.complete,ccc.ward),
			     "Dissimilarity"=dm)
    res <- rbind(res,tmp)
}
kable(res)

#Average Linkage is best independt of dissim measure
clustDat <- filter(clustDat, Statistic=="Average Silhouette Width")
compPlot <- ggplot(clustDat, aes(x=as.character(DeepSplit), y=Value, color=Linkage,
			    group=Linkage)) +
geom_point() +
geom_path() +
facet_grid(.~Dissimilarity, scales="free_y")

compPlot

#Compare Cluster stability 
library(fpc)
boot1 <- bootDat[["euclidean_average_0"]]

#Clustering
clust1 <- dynamicCluster(m.sub, dm="euclidean",
			 lk="average",
			 ds=0,
			 output="Cluster")
clust2 <- dynamicCluster(m.sub, dm="pearson",
			 lk="average",
			 ds=0,
			 output="Cluster")
clust3 <- dynamicCluster(m.sub, dm="pearson",
			 lk="ward.D2",
			 ds=0,
			 output="Cluster")

#Looking a different clusters in depth
pD$cluster <- clust1$cluster
distance <- clust1$dis

ord <- arrange(pD, cluster) 
simMat <- 1-(distance-min(distance))/(max(distance)-min(distance))
simMat.ord <- as.matrix(simMat)[as.character(ord$barcode),
			 as.character(ord$barcode)]
library(pheatmap)
annoCol <- data.frame("cluster"=as.factor(ord$cluster),
		      "Condition"=ord$Condition)
rownames(annoCol) <- as.character(ord$barcode)
pheatmap(simMat.ord,
	 cluster_rows=FALSE,
	 cluster_cols=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE,
	 annotation_col=annoCol)

fPCA <- t(trafM)
fPCA <- scale(fPCA,scale=TRUE,center=TRUE)
set.seed(rnd_seed)
tsn <- Rtsne(fPCA,perplexity=10,initial_dims=5)
pD$tSNE1 <- tsn$Y[,1]
pD$tSNE2 <- tsn$Y[,2]
tsnePlot <- ggplot(pD, aes(x=tSNE1, y=tSNE2, shape=Condition, color=as.factor(cluster))) +
    geom_point(alpha=0.8,size=2) +
    theme_bw()
tsnePlot

tsnePlot <- ggplot(pD.filtered, aes(x=tSNE1, y=tSNE2, color=log2(Acta2+1), shape=Replicate)) +
    geom_point(alpha=0.8,size=2) +
    theme_bw()
tsnePlot

cv2 <- function(x) {
    # Computes cv2 for a vector
    out <- (sd(x)/mean(x))^2
}

normAndDimR <- function(m, pD, fD, top=0.1, scale_true=TRUE) {

	#Thresholding on lowly expressed genes
	m <- m[fD$keep,]
	fD <- fD[fD$keep,]

	clusters <- quickCluster(m)
	pD$sf <- computeSumFactors(m,clusters=clusters,positive=TRUE)
	m.norm <- t(t(m)/pD$sf)
	m.norm <- m.norm[,!is.na(colSums(m.norm))]
	pD <- filter(pD, as.character(barcode) %in% colnames(m.norm))
	# 

	fit <- trendVar(log10(m.norm+1), trend="loess")
	decVar <- decomposeVar(log10(m.norm+1),fit=fit)

	fD <- mutate(fD, meanExpression=rowMeans(m.norm),
		     CV2=apply(m.norm,1,cv2),
		     dm=DM(meanExpression,CV2))

	highVar <- dplyr::arrange(fD,desc(dm)) %>% .$id
	topX <- top*length(highVar)
	fD$highVar <- fD$id %in% highVar[1:topX]
	#         fD$highVar <- fD$id %in% rownames(decVar[decVar$FDR <0.05,])

	fPCA <- log2(t(m.norm[fD$highVar,])+1)
	fPCA <- scale(fPCA,scale=scale_true,center=TRUE)
	pc <- prcomp(fPCA)
	set.seed(rnd_seed)
	tsn <- Rtsne(fPCA,perplexity=floor(ncol(fPCA)/10))

	Acta2 <- filter(fD, symbol=="Acta2") %>% .$id
	pD$Acta2 <- t(m.norm)[,Acta2]
	Cd24 <- filter(fD, symbol=="Cd24a") %>% .$id
	pD$Cd24 <- t(m.norm)[,Cd24]
	Gata3 <- filter(fD, symbol=="Gata3") %>% .$id
	pD$Gata3 <- t(m.norm)[,Gata3]
	Csn2 <- filter(fD, symbol=="Csn2") %>% .$id
	pD$Csn2 <- t(m.norm)[,Csn2]

	pD$PC1 <- pc$x[,1]
	pD$PC2 <- pc$x[,2]
	pD$PC3 <- pc$x[,3]
	pD$tSNE1 <- tsn$Y[,1]
	pD$tSNE2 <- tsn$Y[,2]
	out <- pD
	return(out)
}

plotGeneDist <- function(m, pD, fD, genes, colorBy="Condition") {
    require(reshape)
    stopifnot(identical(rownames(m),fD$id))
    rownames(m) <- fD$symbol
    feat <- filter(fD, id %in% genes) %>% .$symbol
    if (length(feat)==0) {
	feat <- genes
    }
    exps <- data.frame(t(m)[,feat])
    exps$barcode <- rownames(exps)
    pltDat <- left_join(select_(pD, "barcode", colorBy),exps) %>% melt(.,id=c("barcode",colorBy))
    pltDat$value <- log2(pltDat$value+1)
    pltDat[,colorBy] <- as.factor(pltDat[,colorBy])
    plt <- ggplot(pltDat, aes_string(x=colorBy,y="value")) +
	geom_boxplot() +
	geom_point(position="jitter",alpha=0.2,shape=19) +
	facet_wrap(~variable, scales="free") +
	theme_bw()
    return(plt)
}

compClustering <- function(m,fD,fs="brennecke",dm="spearman",lk="ward.D2",ds=0,minSize=15) {

    #Feature Selection
    keep <- fD[,fs]
    trafM <- log2(m[keep,]+1)

    #Dissimilarity Measure
    if(dm!="euclidean") {
	if(lk=="ward.D2") {
    dis <- as.dist(sqrt((1-cor(trafM,method=dm))/2))
	}else {
	    dis <- as.dist((1-cor(trafM,method=dm))/2)
	}
    } else {
	dis <- dist(t(trafM),method=dm)
    }

    #Hierarchical Clustering
    htree <- hclust(dis, method=lk)
    #Tree Cut
    cluss <- cutreeDynamic(htree,distM=as.matrix(dis),deepSplit=ds,minClusterSize=minSize)

    #Compute Validation Statistics
    clstats <- cluster.stats(d=dis, cluss)
    asw <- clstats$avg.silwidth
    rss <- clstats$within.cluster.ss
    con <- connectivity(distance=dis,clusters=cluss)

    #Compile output
    out <- data.frame("FeatureSelection"=fs,
		      "Dissimilarity"=dm,
		      "Linkage"=lk,
		      "DeepSplit"=ds,
		      "Statistic"=c("Average Silhouette Width",
				    "Within Cluster SS",
				    "Connectivity"),
		      "Value"=c(asw,rss,con),
		      "K"=max(cluss),
		      "m"=length(cluss))
}

BrenneckeHVG <- function (m, suppress.plot = FALSE, fdr = 0.1, 
    minBiolDisp = 0.5) 
    {
	require(statmod)
	meansGenes <- rowMeans(m, na.rm = TRUE)
	varsGenes <- unlist(apply(m, 1, var, na.rm = TRUE))
	cv2Genes <- varsGenes/meansGenes^2
	minMeanForFit <- unname(quantile(meansGenes[which(cv2Genes > 0.3)], 
	    0.8))
	useForFit <- meansGenes >= minMeanForFit
	fit <- glmgam.fit(cbind(a0 = 1, a1tilde = 1/meansGenes[useForFit]), 
	    cv2Genes[useForFit])
	a0 <- unname(fit$coefficients["a0"])
	a1 <- unname(fit$coefficients["a1tilde"])
	psia1theta <- a1
	minBiolDisp <- minBiolDisp^2
	m <- ncol(m)
	cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
	testDenom <- (meansGenes * psia1theta + meansGenes^2 * cv2th)/(1 + 
	    cv2th/m)
	p <- 1 - pchisq(varsGenes * (m - 1)/testDenom, m - 1)
	padj <- p.adjust(p, "BH")
	sig <- padj < fdr
	sig[is.na(sig)] <- FALSE
	if (!suppress.plot) {
	    plot(meansGenes, cv2Genes, xaxt = "n", yaxt = "n", log = "xy", 
		xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)", 
		col = "white")
	    axis(1, 10^(-2:5), c("0.01", "0.1", "1", "10", "100", 
		"1000", expression(10^4), expression(10^5)))
	    axis(2, 10^(-2:3), c("0.01", "0.1", "1", "10", "100", 
		"1000"), las = 2)
	    abline(h = 10^(-2:1), v = 10^(-1:5), col = "#D0D0D0", 
		lwd = 2)
	    points(meansGenes, cv2Genes, pch = 20, cex = 0.2, col = ifelse(padj < 
		0.1, "#C0007090", "black"))
	    xg <- 10^seq(-2, 6, length.out = 1000)
	    lines(xg, (a1)/xg + a0, col = "#FF000080", lwd = 3)
	    lines(xg, psia1theta/xg + a0 + minBiolDisp, lty = "dashed", 
		col = "#C0007090", lwd = 3)
	}
	return(names(meansGenes)[sig])
    }
dynamicCluster <- function(m, dm="euclidean", lk="average", ds=0, output="ForBootstrap", minSize=15) {
    print("Clustering on n*p matrix")
    trafM <- log2(m+1)

    #Dissimilarity Measure
    if(dm!="euclidean") {
	if(lk=="ward.D2") {
    dis <- as.dist(sqrt((1-cor(t(trafM),method=dm))/2))
	}else {
	    dis <- as.dist((1-cor(t(trafM),method=dm))/2)
	}
    } else {
	dis <- dist(trafM,method=dm)
    }


    #Hierarchical Clustering
    htree <- hclust(dis, method=lk)
    #Tree Cut
    cluss <- cutreeDynamic(htree,distM=as.matrix(dis),deepSplit=ds,minClusterSize=minSize)

    if (output=="ForBootstrap") {
	out <- list()
	#Sort cluster numbers so that the noise component (0) is always the last cluster (required for cllist)

	csids <- sort(unique(unname(cluss)),decreasing=FALSE)
	out$nc <- length(csids)

	out$clusterlist <- lapply(csids, function(id) cluss %in% id)
	#         names(out$clusterlist) <- csids
	out$partition <- unname(cluss)
	out$clustermethod <- "DynamicTreeCut"
	return(out)
    } else {
	out <- list()
	out$cluster <- cluss
	out$tree <- htree
	out$dis <- dis
	return(out)
    }
}

htreeWithBars <- function(obj,pD)
{
    #takes an obj output from dynamicCluster
    require(dendextend)
    dend <- obj$tree %>% as.dendrogram %>% set("hang")
    clusters <- obj$cluster[order.dendrogram(dend)]
    clusters_numbers <- unique(clusters) - (0 %in% clusters)
    n_clusters <- length(clusters_numbers)
    clusters <- factor(clusters)
    plot(dend,leaflab="none")
    condition <- pD$Condition[order.dendrogram(dend)]
    rplicate <- pD$Replicate[order.dendrogram(dend)]
    levels(condition) <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
    colorBar <- data.frame(clusters,condition,rplicate)
    colored_bars(colorBar, dend, sort_by_labels_order = FALSE, y_shift=-max(obj$tree$height)*0.05,
	     y_scale=max(obj$tree$height)*0.1,rowLabels=c("Cluster","Condition","Animal"))
}

asSim <- function(distance) {
    simMat <- 1-(distance-min(distance))/(max(distance)-min(distance))
}

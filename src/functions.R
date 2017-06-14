# functions used in the analysis
plotGeneDist <- function(m, pD, fD, genes, colorBy="Condition") {
    # function to plot gene expression as boxplots
    require(reshape)
    stopifnot(identical(rownames(m),fD$id))
    rownames(m) <- fD$symbol
    feat <- filter(fD, id %in% genes) %>% .$symbol
    if (length(feat)==0) {
	feat <- genes
    }
    exps <- data.frame(t(m)[,feat])
    colnames(exps) <- feat
    exps$barcode <- rownames(exps)
    pltDat <- left_join(select_(pD, "barcode", colorBy),exps) %>% melt(.,id=c("barcode",colorBy))
    pltDat$value <- log2(pltDat$value+1)
    pltDat[,colorBy] <- as.factor(pltDat[,colorBy])
    plt <- ggplot(pltDat, aes_string(x=colorBy,y="value")) +
	geom_boxplot() +
	#         geom_point(position="jitter",alpha=0.2,shape=19) +
	facet_wrap(~variable, scales="free") +
	theme_bw()
    return(plt)
}

BrenneckeHVG <- function (m, suppress.plot = FALSE, fdr = 0.1, 
    minBiolDisp = 0.25) 
    {
	require(statmod)
	means <- rowMeans(m, na.rm = TRUE)
	vars <- unlist(apply(m, 1, var, na.rm = TRUE))
	cv2 <- vars/means^2
	minMeanForFit <- unname(quantile(means[which(cv2 > 0.3)], 
	    0.8))
	useForFit <- means >= minMeanForFit
	fit <- glmgam.fit(cbind(a0 = 1, a1tilde = 1/means[useForFit]), 
	    cv2[useForFit])
	a0 <- unname(fit$coefficients["a0"])
	a1 <- unname(fit$coefficients["a1tilde"])
	psia1theta <- a1
	minBiolDisp <- minBiolDisp^2
	m <- ncol(m)
	cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
	testDenom <- (means * psia1theta + means^2 * cv2th)/(1 + 
	    cv2th/m)
	p <- 1 - pchisq(vars * (m - 1)/testDenom, m - 1)
	padj <- p.adjust(p, "BH")
	sig <- padj < fdr
	sig[is.na(sig)] <- FALSE
	if (!suppress.plot) {
	    plot(means, cv2, xaxt = "n", yaxt = "n", log = "xy", 
		xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)", 
		col = "white")
	    axis(1, 10^(-2:5), c("0.01", "0.1", "1", "10", "100", 
		"1000", expression(10^4), expression(10^5)))
	    axis(2, 10^(-2:3), c("0.01", "0.1", "1", "10", "100", 
		"1000"), las = 2)
	    abline(h = 10^(-2:1), v = 10^(-1:5), col = "#D0D0D0", 
		lwd = 2)
	    points(means, cv2, pch = 20, cex = 0.2, col = ifelse(padj < 
		0.1, "#C0007090", "black"))
	    xg <- 10^seq(-2, 6, length.out = 1000)
	    lines(xg, (a1)/xg + a0, col = "#FF000080", lwd = 3)
	    lines(xg, psia1theta/xg + a0 + minBiolDisp, lty = "dashed", 
		col = "#C0007090", lwd = 3)
	}
	return(names(means)[sig])
    }
dynamicCluster <- function(dis, lk="average", ds=0, output="ForBootstrap", minSize=15) {
    # ensure dis is dis matrix
    dis <- as.dist(dis)
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

asDis <- function(m, dm="euclidean") {
		logM <- t(log2(m+1))
		#Dissimilarity Measure
		if(dm=="euclidean") {
		    dis <- dist(logM,method=dm)
		} else {
		    dis <- as.dist((1-cor(t(logM),method=dm))/2)
		}
		return(dis)
}


grab_grob <- function(){
  grid.echo()
  grid.grab()
}
compClustering <- function(m,fD,fs="brennecke",dm="spearman",lk="ward.D2",ds=0,minSize=15) {

    #Feature Selection
    keep <- fD[,fs]
    trafM <- t(log2(m[keep,]+1))

    #Dissimilarity Measure
    if(dm=="euclidean") {
	dis <- dist(trafM,method=dm)
    } else {
	dis <- as.dist((1-cor(t(trafM),method=dm))/2)
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


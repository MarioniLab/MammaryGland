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

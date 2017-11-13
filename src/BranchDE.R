# Differential expression along the branches

library(plyr)
library(splines)
library(dplyr)
source("functions.R")
library(lmtest)
library(MASS)

dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

dms <- read.csv("../data/Robjects/dm_luminal.csv")
pD <- right_join(pD,dms,by="barcode")

#sbst m
m <- m[,as.character(pD$barcode)]

#Genes
keep <- rowMeans(m) > 0.1
m <- m[keep,]
fD <- fD[keep,]

# Normalize
m.norm <- t(t(m)/pD$sf)
# ---- BranchDE ----

lineages <- c("Hormone-sensing lineage",
	      "Secretory lineage")
out <- list()
for (lin in lineages) {
    #Subset to lineage+root
    pD.sub <- pD[pD$branch %in% c("Root","Intermediate", 
			      lin),]
    pD.sub$DPTRank <- rank(pD.sub$dpt, ties.method="first")
    m.sub <- m.norm[,as.character(pD.sub$barcode)]
    keep <- rowMeans(m.sub)>0.1 # smoothed values are computed for all, genes are removed afterwards
    m.sub <- log2(m.sub+1)
    fD.sub <- fD

    ids <- colnames(pD.sub)
    yhet <- data.frame(t(m.sub))

    #change ensemblIDs to geneSymbol
    stopifnot(identical(rownames(m.sub),as.character(fD.sub$id)))
    rownames(m.sub) <- colnames(yhet) <- genes <- fD.sub$symbol
    yhet$barcode <- colnames(m.sub)
    fullDat <- join(pD.sub,yhet, by="barcode")

    #initialize m.smooth matrix
    m.smooth <- matrix(nrow=length(genes),ncol=nrow(fullDat))
    colnames(m.smooth) <- as.character(fullDat$barcode)
    rownames(m.smooth) <- genes
    # take m.smooth het to allow allocation of vectors for columns
    m.smooth <- t(m.smooth)

    #initialize result df
    res <- data.frame()
    for (gene in genes) {
	# Set x and y
	x <- fullDat$dpt
	y <- fullDat[,gene]

	# Define null and alternative model
	mod0 <- lm(y ~ 1)
	mod1 <- lm(y ~ ns(x,df=3))

	# Extract coefficients
	cfs <- mod1$coefficients
	names(cfs) <- paste0("c",c(0:(length(cfs)-1)))

	# Likelihood ratio test
	lrt <- lrtest(mod0,mod1)
	p <- lrt[2,5]

	# Linear model for gradient
	lmmod <- lm(mod1$fitted.values~x)
	pgrad <- summary(lmmod)[[4]][2,4]
	gradient <- ifelse(pgrad < 0.01, lmmod$coefficients[2],0)

	# Combine in df
	tmp <- data.frame(Gene=gene,
			  PValue=p,
			  gradient=gradient)
	tmp <- cbind(tmp,t(cfs))
	res <- rbind(res,tmp)

	# Update fitted Value gene expression matrix
	m.smooth[,gene] <- mod1$fitted.values
    }

    # m.smooth back to p*n
    m.smooth <- t(m.smooth)

    # Adjust for multiple testing
    keep <- keep & ((apply(m.smooth,1,max)-apply(m.smooth,1,min))>=0.5)
    res <- res[keep,]
    res$PAdjust<- p.adjust(res$PValue)

    # Store all results in one list per branch
    ord <- arrange(pD.sub, DPTRank) %>% .$barcode %>% as.character()
    m.smooth <- m.smooth[,ord]
    out[[lin]] <- list("Results"=res,
		       "mSmooth"=m.smooth,
		       "m"=m.sub,
		       "pD"=pD.sub)
}

#Combine results from both branches in one DF
hrm <- out[[1]][["Results"]]
colnames(hrm) <- c("Gene",paste0("hrm.",colnames(hrm)[-1]))

alv <- out[[2]][["Results"]]
colnames(alv) <- c("Gene",paste0("alv.",colnames(alv)[-1]))
res <- full_join(hrm,alv,id="Gene")

## Set gradients to 0 if NA
res[,c("alv.gradient","hrm.gradient")][is.na(res[,c("alv.gradient","hrm.gradient")])] <- 0

## Set Pvalues to 1 if NA
res[,c("alv.PAdjust","hrm.PAdjust")][is.na(res[,c("alv.PAdjust","hrm.PAdjust")])] <- 1

# ---- BranchSpecificDefinition ----

#
# Set1 DE on both same gradient
#

# DE definition
res.sameGrad <- filter(res, (hrm.PAdjust < 0.01 & alv.PAdjust < 0.01)
	       & (sign(hrm.gradient)==sign(alv.gradient)))
res.sameGrad <- mutate(res.sameGrad, pValRank=rank(pmin(hrm.PAdjust,alv.PAdjust),ties.method="first"))
out[["genes.sameGrad"]] <- arrange(res.sameGrad, pValRank) %>% .$Gene


#
# Set2 DE with different trends
#

res.diffGrad <- filter(res, (hrm.PAdjust < 0.01 | alv.PAdjust < 0.01)
	       & (sign(hrm.gradient)!=sign(alv.gradient)))

res.diffGrad <- mutate(res.diffGrad, pValRank=rank(pmin(hrm.PAdjust,alv.PAdjust),ties.method="first"))
out[["genes.diffGrad"]] <- arrange(res.diffGrad, pValRank) %>% .$Gene

# Save DE tables

table.sameGrad <-res.sameGrad[,!grepl("c",colnames(res.sameGrad))]
table.sameGrad$hrm.gradient <- sign(table.sameGrad$hrm.gradient)
table.sameGrad$alv.gradient <- sign(table.sameGrad$alv.gradient)
colnames(table.sameGrad) <- gsub("hrm","HormoneSensing",colnames(table.sameGrad))
colnames(table.sameGrad) <- gsub("alv","Secretory",colnames(table.sameGrad))
write.csv(table.sameGrad,"../data/DE_sameGradient.csv",quote=FALSE,row.names=FALSE)

table.diffGrad <-res.diffGrad[,!grepl("c",colnames(res.diffGrad))]
table.diffGrad$hrm.gradient <- sign(table.diffGrad$hrm.gradient)
table.diffGrad$alv.gradient <- sign(table.diffGrad$alv.gradient)
colnames(table.diffGrad) <- gsub("hrm","HormoneSensing",colnames(table.diffGrad))
colnames(table.diffGrad) <- gsub("alv","Secretory",colnames(table.diffGrad))
write.csv(table.diffGrad,"../data/DE_diffGradient.csv",quote=FALSE,row.names=FALSE)

saveRDS(out,"../data/Robjects/BranchDEList.rds")

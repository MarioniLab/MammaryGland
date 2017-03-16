library(plyr)
library(pheatmap)
library(reshape2)
library(splines)
library(dplyr)
library(destiny)
library(ggplot2)
source("functions.R")
library(RColorBrewer)
library(lmtest)

#REad in data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#load transcription factors
tfCheck <- read.table("../data/miscData/TFcheckpoint_WithENSID.tsv",
		header=TRUE, sep="\t")

tfs <- filter(fD, id %in% tfCheck$ensembl_gene_id) %>% .$symbol

# Pre-Filtering  selecting only virgin and pregnancy
# Cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier & !(pD$cluster %in% c(6,7,10)) &
    pD$Condition %in% c("V","P")
m <- m[,keepCells]
pD <- pD[keepCells,]

#Gene filtering
keep <- rowMeans(m)>0.1
m <- m[keep,]
fD <- fD[keep,]

#Normalize
m.norm <- t(t(m)/pD$sf)

# Selection of HVGs
brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
fD$highVar <- fD$id %in% brennecke

#Expression matrix for DiffusionMap
exps <- m.norm[fD$highVar,]
exps <- t(log(exps+1))

# Compute Diffusion map
set.seed(rnd_seed)
dm <- DiffusionMap(exps,n_eigs=20,k=50)
plot(eigenvalues(dm)[1:20])

# Pseudotime and branching
set.seed(rnd_seed)
dcs <- eigenvectors(dm)[,1:2]
t1 <- which.min(dcs[,2])
t2 <- which.min(dcs[,1])
t3 <- which.max(dcs[,1])
dpt <- DPT(dm, branching=TRUE, tips=c(t1,t2,t3))
root <- which(dpt@tips[,1])[3]
rootdpt <- paste0("DPT",root)
plot(dpt,pch=20,dcs=c(1,2),col_by=rootdpt)

# Rename branches
branch <- dpt@branch[,1]
branch[is.na(branch)]  <- "Decision stage"
branch[branch==1] <- "Root"
branch[branch==2] <- "Secretory lineage"
branch[branch==3] <- "Hormone-sensing lineage"
pD$branch <- factor(branch,levels=c("Root","Decision stage",
				    "Secretory lineage",
				    "Hormone-sensing lineage"))
pD$dpt<- dpt[["dpt"]]
pD$DC1 <- eigenvectors(dm)[,1]
pD$DC2 <- eigenvectors(dm)[,2]


#Sanity checks
#PDT is correlated with real time
boxp <- ggplot(pD, aes(SampleID, dpt,fill=Condition)) +
    geom_boxplot() +
    theme_bw(base_size=14)

# Branches
branchPlot <- ggplot(pD, aes(DC1,DC2, color=branch)) +
    geom_point() +
    theme_bw(base_size=14)


#Branch interpretation follwos dpt
boxp2 <- ggplot(pD, aes(branch, dpt)) +
    geom_boxplot() +
    theme_bw(base_size=14)

########################
#
# Hormone-sensing lineage
#
########################

lin <- "Hormone-sensing lineage"
pD.sub <- pD[pD$branch %in% c("Root","Decision stage", 
			      lin),]
pD.sub$DPTRank <- rank(pD.sub$dpt, ties.method="first")
m.sub <- m.norm[,as.character(pD.sub$barcode)]
keep <- rowMeans(m.sub) > 0.01
m.sub <- m.sub[keep,]
fD.sub <- fD[keep,]
ids <- colnames(pD.sub)
yhet <- data.frame(t(m.sub))
stopifnot(identical(rownames(m.sub),as.character(fD.sub$id)))
colnames(yhet) <- fD.sub$symbol
yhet$barcode <- colnames(m.sub)
fullDat <- join(pD.sub,yhet, by="barcode")

res <- data.frame()
genes <- fD.sub$symbol
m.smooth <- matrix(nrow=length(genes),ncol=nrow(fullDat))
colnames(m.smooth) <- as.character(fullDat$barcode)
rownames(m.smooth) <- genes
# take m.smooth het to allow allocation of vectors for columns
m.smooth <- t(m.smooth)
for (gene in genes) {
    x <- fullDat$DPTRank
    y <- log2(fullDat[,gene]+1)
    mod1 <- lm(y ~ ns(x,df=5))
    lmmod <- lm(y~x)

    # F-test only on alternative model
    fstat <- summary(mod1)$fstatistic
    cfs <- mod1$coefficients
    names(cfs) <- paste0("c",c(0:(length(cfs)-1)))
    ## Create P-Value and coefficient matrix
    p <- pf(fstat[1],df1=fstat[2],df2=fstat[3], lower.tail=FALSE)
    gradient <- lmmod$coefficients[2]
    
    tmp <- data.frame(Gene=gene,
		      PValue=p,
		      gradient=gradient)
    tmp <- cbind(tmp,t(cfs))
    res <- rbind(res,tmp)
    ##Update fitted Value gene expression matrix
    m.smooth[,gene] <- mod1$fitted.values
}
m.smooth <- t(m.smooth)

res$PAdjust<- p.adjust(res$PValue, method="bonferroni")

deGenes <- as.character(res$Gene[res$PAdjust< 0.001 &
			res$gradient > 0.0005])
# deGenes <- deGenes[deGenes %in% tfs]


#Prepare heatmap
ord <- arrange(pD.sub, DPTRank) %>% .$barcode %>% as.character()
m.smooth.ord <- m.smooth[deGenes,ord]
rownames(m.sub) <- fD.sub$symbol
m.sub.ord <- log2(m.sub[deGenes,ord]+1)
m.smooth.ord <- m.smooth.ord/apply(m.smooth.ord,1,max)

pheatmap(m.smooth.ord,
	 cluster_cols=FALSE,
	 show_colnames=FALSE,
	 show_rownames=FALSE)

features <- c("Bcl11a","Krt8","Aldh1a3",
	      "Cd14")
forPlot <- fullDat[,c("DPTRank","dpt","barcode","cluster",features)]
forPlot <- melt(forPlot,id=c("barcode","DPTRank","dpt","cluster"))
ggplot(forPlot, aes(x=DPTRank, y=log2(value+1))) +
    geom_point() +
    geom_smooth(method="lm",formula="y~ns(x,df=5)",aes(color=NULL)) +
    facet_wrap(~variable,scales="free") +
    theme_bw()

########################
#
# Two Branch comparison
#
########################

pD.sub <- pD[(!pD$branch %in% c("Root","Decision stage")),] 
pD.sub <- group_by(pD.sub, branch) %>%
    mutate(dptNorm=dpt/max(dpt)) %>%
    ungroup()
			      
m.sub <- m.norm[,as.character(pD.sub$barcode)]
keep <- rowMeans(m.sub) > 0.01
m.sub <- m.sub[keep,]
fD.sub <- fD[keep,]
ids <- colnames(pD.sub)
yhet <- data.frame(t(m.sub))
stopifnot(identical(rownames(m.sub),as.character(fD.sub$id)))
colnames(yhet) <- fD.sub$symbol
yhet$barcode <- colnames(m.sub)
fullDat <- join(pD.sub,yhet, by="barcode")

res <- data.frame()
genes <- fD.sub$symbol
m.smooth <- matrix(nrow=length(genes),ncol=nrow(fullDat))
colnames(m.smooth) <- as.character(fullDat$barcode)
rownames(m.smooth) <- genes
# take m.smooth het to allow allocation of vectors for columns
m.smooth <- t(m.smooth)
for (gene in genes) {
    x <- fullDat$dptNorm
    y <- log2(fullDat[,gene]+1)
    fctr <- fullDat$branch
    mod0 <- lm(y ~ ns(x,df=5))
    mod1 <- lm(y ~ ns(x,df=5)*fctr)
    lrt <- lrtest(mod0,mod1)

    ## Create P-Value and coefficient matrix
    p <- lrt[5][2,]
    
    tmp <- data.frame(Gene=gene,
		      PValue=p)
    res <- rbind(res,tmp)
    ##Update fitted Value gene expression matrix
    #     m.smooth[,gene] <- mod1$fitted.values
}
# m.smooth <- t(m.smooth)

res$PAdjust<- p.adjust(res$PValue, method="bonferroni")

deGenes <- as.character(res$Gene[res$PAdjust< 0.001])
deGenes <- arrange(res, PValue) %>% .$Gene %>% .[1:10]
# deGenes <- deGenes[deGenes %in% tfs]


#Prepare heatmap
ord <- arrange(pD.sub, DPTRank) %>% .$barcode %>% as.character()
m.smooth.ord <- m.smooth[deGenes,ord]
rownames(m.sub) <- fD.sub$symbol
m.sub.ord <- log2(m.sub[deGenes,ord]+1)
m.smooth.ord <- m.smooth.ord/apply(m.smooth.ord,1,max)

pheatmap(m.smooth.ord,
	 cluster_cols=FALSE,
	 show_colnames=FALSE,
	 show_rownames=FALSE)

features <- as.character(deGenes)
forPlot <- fullDat[,c("dptNorm","barcode","branch",features)]
forPlot <- melt(forPlot,id=c("barcode","dptNorm","branch"))
ggplot(forPlot, aes(x=dptNorm, y=log2(value+1),lty=branch)) +
    geom_point() +
    geom_smooth(method="lm",formula="y~ns(x,df=5)",aes(color=NULL)) +
    facet_wrap(~variable,scales="free") +
    theme_bw()

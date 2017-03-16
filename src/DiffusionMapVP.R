library(plyr)
library(reshape2)
library(splines)
library(dplyr)
library(destiny)
library(ggplot2)
source("functions.R")
library(RColorBrewer)
library(e1071)
library(dynamicTreeCut)
library(rgl)

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

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

exps <- m.norm[fD$highVar,]
exps <- t(log(exps+1))

# Compute Diffusion map)
set.seed(rnd_seed)
dm <- DiffusionMap(exps,n_eigs=20,k=50)
plot(eigenvalues(dm)[1:20])

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

#Analyse hormone-sensing lineage
lin <- "Hormone-sensing lineage"
pD.sub <- pD[pD$branch %in% c("Root","Decision stage", 
			      lin),]
pD.sub$DPTRank <- rank(pD.sub$dpt, ties.method="first")
m.sub <- m.norm[,as.character(pD.sub$barcode)]
keep <- rowSums(m.sub>0) > 10
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
    mod <- lm(y ~ ns(x,df=5))
    fstat <- summary(mod)$fstatistic
    cfs <- mod$coefficients
    names(cfs) <- paste0("c",c(0:(length(cfs)-1)))
    ## Create P-Value and coefficient matrix
    p <- pf(fstat[1],df1=fstat[2],df2=fstat[3], lower.tail=FALSE)
    tmp <- data.frame(Gene=gene,
		      PValue=p)
    tmp <- cbind(tmp,t(cfs))
    res <- rbind(res,tmp)
    ##Update fitted Value gene expression matrix
    m.smooth[,gene] <- mod$fitted.values
}
m.smooth <- t(m.smooth)
res$PAdjust <- p.adjust(res$PValue,method="bonferroni")
deGenes <- as.character(res$Gene[res$PAdjust < 0.001])
deGenes <- as.character(arrange(res, PValue) %>% .$Gene %>% .[1:100])




ord <- arrange(pD.sub, DPTRank) %>% .$barcode %>% as.character()
m.smooth.ord <- m.smooth[deGenes,ord]
m.smooth.ord <- m.smooth.ord/apply(m.smooth.ord,1,max)

library(pheatmap)
cairo_pdf("test.pdf",width=12,height=12)
pheatmap(m.smooth.ord,
	 cluster_cols=FALSE,
	 show_colnames=FALSE,
	 show_rownames=TRUE)
dev.off()

rownames(deGenes) <- deGenes$Gene
test <- deGenes[,grepl("^c",colnames(deGenes))]
dis <- (1-cor(t(test)))/2
hcls <- hclust(as.dist(dis))
result <- cutreeDynamic(hcls,distM=as.matrix(dis),deepSplit=0)
names(result) <- deGenes$Gene
# result <- cmeans(test,centers=2)

#cluster genes
rownames(m.sub) <- fD.sub$symbol
m.ord <- m.sub[as.character(deGenes$Gene),pD.sub$DPTRank]
dis <- (1-cor(t(m.ord)))/2
hcls <- hclust(as.dist(dis))
result <- cutreeDynamic(hcls,distM=as.matrix(dis),deepSplit=0)
result <- kmeans(m.ord,2)
result <- result$cluster
names(result) <- deGenes$Gene


# Plot to explore expression of a few genes
c1 <- sample(names(result[result==1]),10)
c2 <- sample(names(result[result==1]),10)
features <- c("Krt8")
forPlot <- fullDat[,c("DPTRank","barcode",features)]
forPlot <- melt(forPlot,id=c("barcode","DPTRank"))
ggplot(forPlot, aes(x=DPTRank, y=log2(value+1))) +
    geom_point() +
    geom_smooth(method="lm",formula="y~ns(x,df=5)") +
    facet_wrap(~variable) +
    theme_bw()



g <- filter(fD, symbol=="") %>% .$id
y <- log((t(m.norm)[,gn])+1)
dat <- data.frame(pseudot,y,cluster=pD$cluster,branch=factor(branch),
		  Condition=pD$Condition) %>%
    filter(!is.na(branch) & branch!=1) %>%
    group_by(branch) %>%
    mutate(DPTRank=rank(pseudot, ties.method="first")) %>%
    ungroup()
ggplot(dat, aes(DPTRank,y,col=cluster)) +
    geom_point() +
    geom_smooth(method="loess",aes(lty=branch),col="black") +
    theme_bw()

ggplot(dat, aes(pseudot, ..density..)) +
    geom_histogram(bins=100) +
    facet_wrap(~branch) +
    theme_bw()

library(pheatmap)
mat <- m.norm[,branch==3 & !is.na(branch)]
keep <- rowMeans(mat) > 0.1
mat <- mat[keep,]
brennecke <- BrenneckeHVG(mat,suppress.plot=TRUE,fdr=.01)
mat <- log2(mat[brennecke,order(pseudot[branch==3 & !is.na(branch)])]+1)
mat <- mat - rowMeans(mat)
rownames(mat) <- fD$symbol[fD$id %in% brennecke]
pheatmap(mat,
	 cluster_cols=FALSE,
	 #          clustering_distance_rows="correlation",
	 #          clustering_method="ward.D2",
	 show_rownames=TRUE,
	 fontsize=6,
	 show_colnames=FALSE)

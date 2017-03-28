library(plyr)
library(cowplot)
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

#Rename clusters
pD$cluster <- mapvalues(pD$cluster, from=c(1,2,3,4,5,6,7,9,10),
			  to=c(1,2,3,4,5,6,7,8,9))

# Pre-Filtering  selecting only virgin and pregnancy
# Cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier & !(pD$cluster %in% c(6,7,9)) &
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

library(viridis)
b1 <- pD[pD$branch %in% c("Root","Decision stage",
			  "Hormone-sensing lineage"),]
pb1 <- ggplot(pD,aes(x=DC1,y=-DC2)) +
    geom_point(size=0.8,color="grey80") +
    geom_point(data=b1,aes(x=DC1,y=-DC2,color=dpt)) +
    scale_color_viridis(option="magma",begin=1,end=0) +
    xlab("Component 1") +
    ylab("-Component 2") +
    ggtitle("Hormone-sensing lineage") +
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank(),
	  legend.title=element_blank())


b2 <- pD[pD$branch %in% c("Root","Decision stage",
			  "Secretory lineage"),]
pb2 <- ggplot(pD,aes(x=DC1,y=-DC2)) +
    geom_point(size=0.8,color="grey80") +
    geom_point(data=b2,aes(x=DC1,y=-DC2,color=dpt)) +
    scale_color_viridis(option="magma",begin=1,end=0) +
    xlab("Component 1") +
    ylab("-Component 2") +
    ggtitle("Secretory lineage") +
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank(),
	  legend.title=element_blank())

branches <- plot_grid(pb2,NULL,pb1,NULL,nrow=1,
		      rel_widths=c(1,.2,1,.4))

########################
#
# Lineage-specific Differential Expression
#
########################
lineages <- c("Hormone-sensing lineage",
	      "Secretory lineage")
out <- list()
for (lin in lineages) {
    pD.sub <- pD[pD$branch %in% c("Root","Decision stage", 
			      lin),]
    pD.sub$DPTRank <- rank(pD.sub$dpt, ties.method="first")
    m.sub <- log2(m.norm[,as.character(pD.sub$barcode)]+1)
    #     keep <- rowSums(m.sub) != 0
    #     m.sub <- m.sub[keep,]
    fD.sub <- fD#[keep,]
    ids <- colnames(pD.sub)
    yhet <- data.frame(t(m.sub))
    stopifnot(identical(rownames(m.sub),as.character(fD.sub$id)))
    rownames(m.sub) <- fD.sub$symbol
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
	x <- fullDat$dpt
	y <- fullDat[,gene]
	mod1 <- lm(y ~ ns(x,df=3))
	lmmod <- lm(mod1$fitted.values~x)

	# F-test only on alternative model
	fstat <- summary(mod1)$fstatistic
	cfs <- mod1$coefficients
	names(cfs) <- paste0("c",c(0:(length(cfs)-1)))
	## Create P-Value and coefficient matrix
	p <- pf(fstat[1],df1=fstat[2],df2=fstat[3], lower.tail=FALSE)
	pgrad <- summary(lmmod)[[4]][2,4]
	gradient <- ifelse(pgrad < 0.01, lmmod$coefficients[2],0)

	tmp <- data.frame(Gene=gene,
			  PValue=p,
			  gradient=gradient)
	tmp <- cbind(tmp,t(cfs))
	res <- rbind(res,tmp)
	##Update fitted Value gene expression matrix
	m.smooth[,gene] <- mod1$fitted.values
    }
    m.smooth <- t(m.smooth)
    maxExp <-apply(m.smooth,1,max)
    minExp <-apply(m.smooth,1,min)
    tooLow <- rownames(m.smooth)[(maxExp-minExp)<=1]


    res$PAdjust<- p.adjust(res$PValue, method="bonferroni")
    res <- mutate(res,tooLow= Gene %in% tooLow)

    ## Heatmap of all TF 
    #     deGenes <- as.character(res$Gene[res$PAdjust< 0.001 & !res$tooLow])
    #     deGenes <- deGenes[deGenes %in% tfs]
    #     deGenes <- as.character(res$Gene[order(res$PValue)][1:50])

    ord <- arrange(pD.sub, DPTRank) %>% .$barcode %>% as.character()
    #     m.heat <- m.smooth[deGenes,ord]
    m.smooth <- m.smooth[,ord]
    #     m.heat <- m.heat/apply(m.heat,1,max)
    #     ords1 <- apply(m.heat,1,function(x) min(which(x>0.5)))
    #     ords2 <- apply(m.heat,1,function(x) min(which(x<0.5)))
    #     ord <- ifelse(ords1==1,ords2,ords1)
    #     ord <- names(sort(ord))
    #     m.heat <- m.heat[ord,]
    #     row.dist <- as.dist((1-cor(t(m.heat),method="spearman"))/2)

    #     pheat <- pheatmap(m.heat,
    #          cluster_cols=FALSE,
    #          cluster_rows=TRUE,
    #          show_colnames=FALSE,
    #          clustering_distance_rows=row.dist,
    #          treeheight_row=0,
    #          clustering_method="ward.D2",
    #          show_rownames=FALSE,
    #          fontsize=6)

    #     if (lin=="Hormone-sensing lineage") {
    #     features <- c("Esr1","Foxa1","Pgr")
    #     lables <- c("a","b")
    #     } else{
    #     features <- c("Csn2","Glycam1")
    #     lables <- c("c","d")
    #     }
    # 
    #     plList <- list()
    #     for (feature in features) {
    #     forPlot <- fullDat[,c("DPTRank","dpt","barcode","cluster",feature)]
    #     forPlot <- melt(forPlot,id=c("barcode","DPTRank","dpt","cluster"))
    #     p <- ggplot(forPlot, aes(x=dpt, y=log2(value+1))) +
    #         geom_point(size=0.7) +
    #         geom_smooth(method="lm",formula="y~ns(x,df=3)") +
    #         xlab("Pseudotiqe") +
    #         ylab("Log-Expression") +
    #         ggtitle(feature) 
    #     plList[[feature]] <- p
    #     }
    #     expPlot <- plot_grid(plotlist=plList)
    # 
    #     subP1 <- plot_grid(expPlot,pTf[[4]],nrow=1,labels=c(lables[1],lables[2]),
    #                        vjust=0.5,scale=0.95)
    #     plt <- plot_grid(pheat[[4]])
    out[[lin]] <- list("Results"=res,
		       "MatrixSmooth"=m.smooth,
		       "Matrix"=m.sub,
		       "pD"=pD.sub)
		       #                        "Plot"=plt)

}

# dev.off()
# cairo_pdf("Figure3.pdf",width=12.41,height=8.75)
# plot_grid(plotlist=list(out[[1]][[1]],out[[2]][[1]]),nrow=1)
# dev.off()

hrm <- out[[1]][["Results"]]
colnames(hrm) <- c("Gene",paste0("hrm.",colnames(hrm)[-1]))

alv <- out[[2]][["Results"]]
colnames(alv) <- c("Gene",paste0("alv.",colnames(alv)[-1]))

res <- inner_join(hrm,alv,id="Gene")
m.hrm <- out[[1]][["MatrixSmooth"]]
# m.hrm <- m.hrm/apply(m.hrm,1,max)
m.alv <- out[[2]][["MatrixSmooth"]]
# m.alv <- m.alv/apply(m.alv,1,max)
genes.both <- intersect(rownames(m.hrm),rownames(m.alv))

m.both <- cbind(m.alv[genes.both,c(ncol(m.alv):1)],m.hrm[genes.both,])
# m.both <- m.both/apply(m.both,1,max)
m.both <- t(scale(t(m.both)))
m.both[m.both > 3]  <- 3
m.both[m.both < -3]  <- -3
# 
#Only DE on both branches
# res <- filter(res, hrm.PAdjust < 0.001 | alv.PAdjust < 0.001)
# res <- filter(res, hrm.tooLow & alv.tooLow)

########### Three set of Genes #################################

#Set1 DE on both same gradient
res1 <- filter(res, (hrm.PAdjust < 0.01 & alv.PAdjust < 0.01)
	       & (sign(hrm.gradient)==sign(alv.gradient)))
	       #                & (!(hrm.tooLow | alv.tooLow)))

cres <- rbind(res1,res1)
combPvalRank <- order(c(res1$hrm.PAdjust,res1$alv.PAdjust))
genes1 <- unique(cres[combPvalRank,"Gene"]) %>%.[1:50] %>% as.character()
swtch1 <- apply(m.both[genes1,],1,function(x) max(which(x>0.5)))
genes1 <- res1$Gene

res1tfs <- filter(res1, Gene %in% tfs) %>% .$Gene %>% as.character()

#Set2 DE with different trends
res2 <- filter(res, (hrm.PAdjust < 0.01 | alv.PAdjust < 0.01)
	       & (sign(hrm.gradient)!=sign(alv.gradient)))
	       #                & (!(hrm.tooLow & alv.tooLow)))
cres <- rbind(res2,res2)
combPvalRank <- order(c(res2$hrm.PAdjust,res2$alv.PAdjust))
genes2 <- unique(cres[combPvalRank,"Gene"]) %>%.[1:50] %>% as.character()
swtch2 <- apply(m.both[genes2,],1,function(x) max(which(x>0.5)))
genes2 <- res2$Gene

res2tfs <- filter(res2, Gene %in% tfs) %>% .$Gene %>% as.character()

forxls1 <-res1[,!grepl("c|tooLow",colnames(res1))]
colnames(forxls1) <- gsub("hrm","HormoneSensing",colnames(forxls1))
colnames(forxls1) <- gsub("alv","Secretory",colnames(forxls1))
write.csv(forxls1,"../paper/supps/DE_sameGradient.csv",quote=FALSE)

forxls2 <-res2[,!grepl("c|tooLow",colnames(res2))]
colnames(forxls2) <- gsub("hrm","HormoneSensing",colnames(forxls2))
colnames(forxls2) <- gsub("alv","Secretory",colnames(forxls2))
write.csv(forxls2,"../paper/supps/DE_diffGradient.csv",quote=FALSE)

    


#Set3 DE increase only hormone-sensing
# res3 <- filter(res, hrm.PAdjust < 0.01 & !hrm.tooLow & hrm.gradient >0)
# &alv.gradient<=0)
# genes3 <- arrange(res3,hrm.PAdjust) %>% .$Gene %>%.[1:30] %>% as.character()
# swtch3 <- apply(m.both[genes3,],1,function(x) min(which(x>0.5)))
# genes3 <- genes3[order(swtch3)]
# 
# res3tfs <- filter(res3, Gene %in% tfs) %>% .$Gene %>% as.character()
# 

m.heat <- m.both

##Annotation Data
pD.ord <- pD
rownames(pD.ord) <- pD.ord$barcode
pD.ord <- pD.ord[colnames(m.heat),]
annoCol <- data.frame("Cluster"=as.factor(pD.ord$cluster),
		      "Pseudotime"=rank(pD.ord$dpt,ties.method="first")
		      )
#Rename the alveolar cells, so that there are no duplicate names for rows
colnames(m.heat) <- c(paste0("Alv.",colnames(m.heat)[1:ncol(m.alv)]),
		       colnames(m.heat)[(ncol(m.alv)+1):ncol(m.heat)])

rownames(annoCol) <- colnames(m.heat)

library(RColorBrewer)
clustCol <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,8)]
names(clustCol) <- c(1,2,3,5,8)
library(viridis)
dptcols <- viridis(n=nrow(annoCol),,option="magma",begin=1,end=0)
names(dptcols) <- c(1:length(dptcols))
annoColors <- list("Cluster"=clustCol,
		   "Pseudotime"=dptcols)

m.heat1 <- m.heat[genes1,]
m.heat2 <- m.heat[genes2,]
##Plot
png("Heatmap_sameGradient.png",height=1248,width=960)
p0 <- pheatmap(m.heat1,
	 cluster_cols=FALSE,
	 cluster_rows=TRUE,
	 clustering_distance_rows="correlation",
	 annotation=annoCol,
	 clustering_method="average",
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 treeheight_row=0,
	 legend=FALSE,
	 annotation_legend=FALSE,
	 gaps_col=ncol(m.alv),
	 show_rownames=FALSE)
dev.off()

png("Heatmap_diffGradient.png",height=1248,width=960)
p1 <- pheatmap(m.heat2,
	 cluster_cols=FALSE,
	 cluster_rows=TRUE,
	 clustering_distance_rows="correlation",
	 annotation=annoCol,
	 clustering_method="average",
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 treeheight_row=0,
	 annotation_legend=TRUE,
	 gaps_col=ncol(m.alv),
	 show_rownames=FALSE)
dev.off()

png("../paper/figures/S6.png",height=1248,width=1248)
comb <- plot_grid(p0[[4]],NULL,p1[[4]],nrow=1,vjust=0.5,rel_widths=c(1,.1,1))
comb
dev.off()

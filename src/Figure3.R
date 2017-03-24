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
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank())

b2 <- pD[pD$branch %in% c("Root","Decision stage",
			  "Secretory lineage"),]
pb2 <- ggplot(pD,aes(x=DC1,y=-DC2)) +
    geom_point(size=0.8,color="grey80") +
    geom_point(data=b2,aes(x=DC1,y=-DC2,color=dpt)) +
    scale_color_viridis(option="magma",begin=1,end=0) +
    xlab("Component 1") +
    ylab("-Component 2") +
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank())

branches <- plot_grid(pb2,pb1)

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
m.both <- m.both/apply(m.both,1,max)
# m.both <- t(scale(t(m.both)))
# m.both[m.both>3] <- 3
# m.both[m.both<-3] <- -3
# 
#Only DE on both branches
# res <- filter(res, hrm.PAdjust < 0.001 | alv.PAdjust < 0.001)
# res <- filter(res, hrm.tooLow & alv.tooLow)

########### Three set of Genes #################################

#Set1 DE on both same gradient
res1 <- filter(res, (hrm.PAdjust < 0.001 & alv.PAdjust < 0.001)
	       & (sign(hrm.gradient)<0 & sign(alv.gradient) < 0)
	       & (!(hrm.tooLow | alv.tooLow)))

cres <- rbind(res1,res1)
combPvalRank <- order(c(res1$hrm.PAdjust,res1$alv.PAdjust))
genes1 <- cres[combPvalRank,"Gene"] %>%.[1:30] %>% as.character()
swtch1 <- apply(m.both[genes1,],1,function(x) max(which(x>0.5)))
genes1 <- genes1[order(swtch1)]

res1tfs <- filter(res1, Gene %in% tfs) %>% .$Gene %>% as.character()

#Set2 DE increase only alveolar
res2 <- filter(res, alv.PAdjust < 0.001 & !alv.tooLow & alv.gradient >0)
#                & hrm.gradient<=0)
genes2 <- arrange(res2,alv.PAdjust) %>% .$Gene %>%.[1:30] %>% as.character()
swtch2 <- apply(m.both[genes2,],1,function(x) max(which(x>0.5)))
genes2 <- genes2[order(swtch2)]

res2tfs <- filter(res2, Gene %in% tfs) %>% .$Gene %>% as.character()



#Set3 DE increase only hormone-sensing
res3 <- filter(res, hrm.PAdjust < 0.001 & !hrm.tooLow & hrm.gradient >0)
	       #                &alv.gradient<=0)
genes3 <- arrange(res3,hrm.PAdjust) %>% .$Gene %>%.[1:30] %>% as.character()
swtch3 <- apply(m.both[genes3,],1,function(x) min(which(x>0.5)))
genes3 <- genes3[order(swtch3)]

res3tfs <- filter(res3, Gene %in% tfs) %>% .$Gene %>% as.character()


genes <- c(genes1,genes2,genes3)
m.heat <- m.both[genes,]

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

##Plot
p0 <- pheatmap(m.heat,
	 cluster_cols=FALSE,
	 cluster_rows=FALSE,
	 clustering_distance_rows="correlation",
	 annotation=annoCol,
	 clustering_method="average",
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 treeheight_row=0,
	 cutree_rows=2,
	 annotation_legend=FALSE,
	 gaps_col=ncol(m.alv),
	 gaps_row=c(30,60),
	 show_rownames=TRUE,
	 fontsize=6)

########################
#
# Two Branch comparison
#
########################


###Try 2
features <- c("Krt8","Krt18","Gata3","Notch1","Notch2","Stat5a","Stat6",
	     "Bmi1","Prom1","Sox9","Elf5","Hes1","Lmo4","Wnt5a","Lalba",
	     "Krt14","Pgr","Tnfsf11")
features <- c("Creb5","Hey1","Fosl1",
	      "Runx1","Tox2","Bhlhe41",
	      "Elf5","Foxs1","Ehf"
	      )
# features <- res1tfs
p1 <- out[[1]][["pD"]] %>% arrange(DPTRank)
yhet1 <- data.frame(t(m.hrm)[,features])
yhet1$barcode <- as.character(p1$barcode)
raw1 <- data.frame(t(out[[1]][["Matrix"]][features,yhet1$barcode]))
colnames(raw1) <- paste0("raw",colnames(raw1))
raw1$barcode <- yhet1$barcode
fplot1 <- join(p1,yhet1,by="barcode") 
fplot1 <- join(fplot1,raw1,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))
p2 <- out[[2]][["pD"]] %>% arrange(DPTRank)
yhet2 <- data.frame(t(m.alv)[,features])
yhet2$barcode <- as.character(p2$barcode)
raw2 <- data.frame(t(out[[2]][["Matrix"]][features,yhet2$barcode]))
colnames(raw2) <- paste0("raw",colnames(raw2))
raw2$barcode <- yhet2$barcode
fplot2 <- join(p2,yhet2,by="barcode") 
fplot2 <- join(fplot2,raw2,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))
pList <- list()
for (feature in features) {
    pnts <- paste0("raw",feature)
    lns <- feature
    p <- ggplot() +
	geom_point(size=0.8,data=fplot1,aes_string(x="dptNorm",y=pnts, color="cluster")) +
	geom_point(size=0.8,data=fplot2,aes_string(x="dptNorm",y=pnts, color="cluster")) +
	geom_line(data=fplot1,aes_string(x="dptNorm",y=lns),lty="dashed") +
	geom_line(data=fplot2,aes_string(x="dptNorm",y=lns)) +
	ggtitle(feature) +
	ylab("") +
	scale_color_manual(values=clustCol) +
	scale_x_continuous(breaks=c(0,1)) +
	theme(legend.position="bottom",
	      legend.direction="horizontal",
	      legend.title=element_blank()) +
	guides(colour = guide_legend(override.aes = list(size=3))) +
	xlab("") 
    pList[[feature]] <- p
}

	      

# pD.sub <- pD[(!pD$branch %in% c("Root","Decision stage")),] 
# pD.sub <- group_by(pD.sub, branch) %>%
#     mutate(dptNorm=dpt/max(dpt)) %>%
#     ungroup()
#                               
# m.sub <- m.norm[,as.character(pD.sub$barcode)]
# keep <- rowMeans(m.sub) > 0.1
# m.sub <- m.sub[keep,]
# fD.sub <- fD#[keep,]
# ids <- colnames(pD.sub)
# yhet <- data.frame(t(m.sub))
# stopifnot(identical(rownames(m.sub),as.character(fD.sub$id)))
# colnames(yhet) <- fD.sub$symbol
# yhet$barcode <- colnames(m.sub)
# fullDat <- join(pD.sub,yhet, by="barcode")

# res <- data.frame()
# genes <- fD.sub$symbol
# m.smooth <- matrix(nrow=length(genes),ncol=nrow(fullDat))
# colnames(m.smooth) <- as.character(fullDat$barcode)
# rownames(m.smooth) <- genes
# take m.smooth het to allow allocation of vectors for columns
# m.smooth <- t(m.smooth)
# for (gene in genes) {
#     x <- fullDat$dptNorm
#     y <- log2(fullDat[,gene]+1)
#     fctr <- fullDat$branch
#     mod0 <- lm(y ~ ns(x,df=5))
#     mod1 <- lm(y ~ ns(x,df=5)*fctr)
#     lrt <- lrtest(mod0,mod1)
# 
    ## Create P-Value and coefficient matrix
#     p <- lrt[5][2,]
#     
#     tmp <- data.frame(Gene=gene,
#                       PValue=p)
#     res <- rbind(res,tmp)
    ##Update fitted Value gene expression matrix
    #     m.smooth[,gene] <- mod1$fitted.values
# }
# m.smooth <- t(m.smooth)
# 
# res$PAdjust<- p.adjust(res$PValue, method="bonferroni")
# 
# deGenes <- as.character(res$Gene[res$PAdjust< 0.001])
# deGenes <- arrange(res, PValue) %>% .$Gene %>% .[1:20]
#  deGenes <- deGenes[deGenes %in% tfs] [1:50]


#Prepare heatmap
# ord <- arrange(pD.sub, DPTRank) %>% .$barcode %>% as.character()
# m.smooth.ord <- m.smooth[deGenes,ord]
# rownames(m.sub) <- fD.sub$symbol
# m.sub.ord <- log2(m.sub[deGenes,ord]+1)
# m.smooth.ord <- m.smooth.ord/apply(m.smooth.ord,1,max)
# 
# pheatmap(m.smooth.ord,
#          cluster_cols=FALSE,
#          show_colnames=FALSE,
#          show_rownames=FALSE)

# forPlot <- fullDat[,c("dptNorm","barcode","branch",features)]
# forPlot <- melt(forPlot,id=c("barcode","dptNorm","branch"))
# ggplot(forPlot, aes(x=dptNorm, y=log2(value+1),color=branch)) +
#     geom_point(size=0.5) +
#     geom_smooth(method="lm",formula="y~ns(x,df=3)") +
#     facet_wrap(~variable,scales="free") +
#     theme_bw()
 
# plList <- list()
# for (feature in features) {
#     forPlot <- fullDat[,c("dptNorm","branch","barcode","cluster",feature)]
#     forPlot <- melt(forPlot,id=c("barcode","dptNorm","branch","cluster"))
#     p <- ggplot(forPlot, aes(x=dptNorm, y=log2(value+1),color=cluster,lty=branch)) +
#         geom_point(size=0.7) +
#         geom_smooth(method="lm",formula="y~ns(x,df=3)",aes(color=NULL),color="black") +
#         xlab("Pseudotime") +
#         ylab("Log-Expression") +
#         ggtitle(feature) +
#         guides(color=FALSE,lty=FALSE)
#     plList[[feature]] <- p
# }

dummyD <- data.frame(x=rnorm(10),y=rnorm(10),
		     branch=c(rep("Secretory lineage",5),rep("Hormone-sensing lineage",5)))
dummyP2 <- ggplot(dummyD,aes(x,y,lty=branch)) +
    geom_line() +
    theme(legend.position="bottom",
	  #           legend.direction="horizontal",
	  legend.title=element_blank()) +
    scale_linetype_manual(values=c("dashed","solid"))
legs1 <- get_legend(dummyP2)

legs <- get_legend(pList[[1]])
legs <- plot_grid(legs1,legs,nrow=1)
pls <- lapply(pList,function(x){
	      x <- x %+% guides(color=FALSE) 
	      return(x)})

expPlot <- plot_grid(plotlist=pls,ncol=3)
expPlot <- plot_grid(expPlot,legs,rel_heights=c(1,0.05),ncol=1)

htmp <- plot_grid(branches,p0[[4]],ncol=1,rel_heights=c(0.5,1))
fullP <- plot_grid(htmp,expPlot,ncol=2,rel_widths=c(1,1))

# save_plot("Figure3.pdf",fullP,base_height=6.7675)
cairo_pdf("Figure3.pdf",width=16.55,height=13.0575)
fullP
dev.off()
# 

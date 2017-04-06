#########################################
#
#Script to reproduce Figure 3
#
#########################################
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
library(viridis)
library(lmtest)
library(gridExtra)

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# ---- OnlyLuminal ----
condComb <- c("NP","G")
# Cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier & !(pD$cluster %in% c(6,7,9)) & pD$Condition %in% condComb
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes 
keep <- rowMeans(m)>0.1
m <- m[keep,]
fD <- fD[keep,]

# Normalize
m.norm <- t(t(m)/pD$sf)

# HVGs
brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
fD$highVar <- fD$id %in% brennecke

# Prepare expression matrix by selecting only HVGs and log-transformation
exps <- m.norm[fD$highVar,]
exps <- t(log(exps+1))

# Compute Diffusion map
set.seed(rnd_seed)
dm <- DiffusionMap(exps,n_eigs=20,k=50)
dcs <- eigenvectors(dm)[,1:2]

#Define tips as the cells at the corners of the triangluar shape
t1 <- which.min(dcs[,2])
t2 <- which.min(dcs[,1])
t3 <- which.max(dcs[,1])

# Compute Pseudotime and branching
set.seed(rnd_seed)
dpt <- DPT(dm, branching=TRUE, tips=c(t1,t2,t3))
root <- which(dpt@tips[,1])[3]
rootdpt <- paste0("DPT",root)

# Rename branches
branch <- dpt@branch[,1]
branch[is.na(branch)]  <- "Intermediate"
branch[branch==1] <- "Root"
branch[branch==2] <- "Secretory lineage"
branch[branch==3] <- "Hormone-sensing lineage"

#add branches and pseudotime to pD
pD$branch <- factor(branch,levels=c("Root","Intermediate",
				    "Secretory lineage",
				    "Hormone-sensing lineage"))
pD$dpt<- dpt[["dpt"]]
pD$DC1 <- eigenvectors(dm)[,1]
pD$DC2 <- eigenvectors(dm)[,2]

# ---- PlotWithBranchDPT -----

b1 <- pD[pD$branch %in% c("Root","Intermediate",
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


b2 <- pD[pD$branch %in% c("Root","Intermediate",
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

# ---- BranchDE ----

lineages <- c("Hormone-sensing lineage",
	      "Secretory lineage")
out <- list()
for (lin in lineages) {
    #Subset to lineage+root
    pD.sub <- pD[pD$branch %in% c("Root","Intermediate", 
			      lin),]
    pD.sub$DPTRank <- rank(pD.sub$dpt, ties.method="first")
    m.sub <- log2(m.norm[,as.character(pD.sub$barcode)]+1)
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
	#Set x and y
	x <- fullDat$dpt
	y <- fullDat[,gene]

	#Define null and alternative model
	mod0 <- lm(y ~ 1)
	mod1 <- lm(y ~ ns(x,df=3))

	#Extract coefficients
	cfs <- mod1$coefficients
	names(cfs) <- paste0("c",c(0:(length(cfs)-1)))

	#Likelihood ratio test
	lrt <- lrtest(mod0,mod1)
	p <- lrt[2,5]

	#Linear model for gradient
	lmmod <- lm(mod1$fitted.values~x)
	pgrad <- summary(lmmod)[[4]][2,4]
	gradient <- ifelse(pgrad < 0.01, lmmod$coefficients[2],0)

	#Combine in df
	tmp <- data.frame(Gene=gene,
			  PValue=p,
			  gradient=gradient)
	tmp <- cbind(tmp,t(cfs))
	res <- rbind(res,tmp)
	##Update fitted Value gene expression matrix
	m.smooth[,gene] <- mod1$fitted.values
    }

    #m.smooth back to p*n
    m.smooth <- t(m.smooth)

    #Adjust for multiple testing
    res$PAdjust<- p.adjust(res$PValue)

    #Store all results in one list per branch
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
res <- inner_join(hrm,alv,id="Gene")

#Combine smoothed expression values for heatmap
m.hrm <- out[[1]][["mSmooth"]]
m.alv <- out[[2]][["mSmooth"]]
m.both <- cbind(m.alv[,c(ncol(m.alv):1)],m.hrm) # reverse order of alveolar cells for heatmap
#Scale expression values
m.both <- t(scale(t(m.both)))
m.both[m.both>3] <- 3 # cut at 3 for visualization
m.both[m.both<-3] <- -3 # cut at -3 for visualization


# ---- BranchSpecificDefinition ----


#Set1 DE on both same gradient
res1 <- filter(res, (hrm.PAdjust < 0.01 & alv.PAdjust < 0.01)
	       & (sign(hrm.gradient)==sign(alv.gradient)))

res1 <- mutate(res1, pValRank=rank(pmin(hrm.PAdjust,alv.PAdjust),ties.method="first"))
genes1 <- arrange(res1, pValRank) %>% .[1:50,"Gene"] %>% as.character()
swtch1 <- apply(m.both[genes1,],1,function(x) max(which(x>0.5)))
genes1 <- genes1[order(swtch1)]


#Set2 DE with different trends
res2 <- filter(res, (hrm.PAdjust < 0.01 | alv.PAdjust < 0.01)
	       & (sign(hrm.gradient)!=sign(alv.gradient)))

res2 <- mutate(res2, pValRank=rank(pmin(hrm.PAdjust,alv.PAdjust),ties.method="first"))
genes2 <- arrange(res2, pValRank) %>% .[1:50,"Gene"] %>% as.character()
swtch2 <- apply(m.both[genes2,],1,function(x) max(which(x>0.5)))
genes2 <- genes2[order(swtch2)]


#Matrix for heatmap
genes <- c(genes1,genes2)
m.heat <- m.both[genes,]

#Annotation Data
pD.ord <- pD
rownames(pD.ord) <- pD.ord$barcode
pD.ord <- pD.ord[colnames(m.heat),]
annoCol <- data.frame("Cluster"=pD.ord$cluster,
		      "Pseudotime"=rank(pD.ord$dpt,ties.method="first")
		      )

#Rename the alveolar cells, so that there are no duplicate names for rows
colnames(m.heat) <- c(paste0("Alv.",colnames(m.heat)[1:ncol(m.alv)]),
		       colnames(m.heat)[(ncol(m.alv)+1):ncol(m.heat)])
rownames(annoCol) <- colnames(m.heat)

#Set colorscheme for heatmap
clustCol <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,8)]
names(clustCol) <- c(1,2,3,5,8)
dptcols <- viridis(n=nrow(annoCol),,option="magma",begin=1,end=0)
names(dptcols) <- c(1:length(dptcols))
annoColors <- list("Cluster"=clustCol,
		   "Pseudotime"=dptcols)

#Plot heatmap
p0 <- pheatmap(m.heat,
	 cluster_cols=FALSE,
	 cluster_rows=FALSE,
	 annotation=annoCol,
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 annotation_legend=FALSE,
	 gaps_col=ncol(m.alv),
	 gaps_row=c(50),
	 show_rownames=TRUE,
	 fontsize=6)


# ---- PlotTFExamples ----


#Set Genes to plot (TFs that are DE)
features <- c("Creb5","Hey1","Fosl1",
	      "Runx1","Tox2","Bhlhe41",
	      "Elf5","Foxs1","Ehf")

#extract expression for first branch
p1 <- out[[1]][["pD"]] %>% arrange(DPTRank)
yhet1 <- data.frame(t(m.hrm)[,features])
yhet1$barcode <- as.character(p1$barcode)
raw1 <- data.frame(t(out[[1]][["m"]][features,yhet1$barcode]))
colnames(raw1) <- paste0("raw",colnames(raw1))
raw1$barcode <- yhet1$barcode
fplot1 <- join(p1,yhet1,by="barcode") 
fplot1 <- join(fplot1,raw1,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))

#extract expression for second branch
p2 <- out[[2]][["pD"]] %>% arrange(DPTRank)
yhet2 <- data.frame(t(m.alv)[,features])
yhet2$barcode <- as.character(p2$barcode)
raw2 <- data.frame(t(out[[2]][["m"]][features,yhet2$barcode]))
colnames(raw2) <- paste0("raw",colnames(raw2))
raw2$barcode <- yhet2$barcode
fplot2 <- join(p2,yhet2,by="barcode") 
fplot2 <- join(fplot2,raw2,by="barcode") %>%
    mutate(dptNorm=dpt/max(dpt))


#Plot dpt-dependent expression for features
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

# create a dummy legend for dasehd/solid line
dummyD <- data.frame(x=rnorm(10),y=rnorm(10),
		     branch=c(rep("Secretory lineage",5),rep("Hormone-sensing lineage",5)))
dummyP2 <- ggplot(dummyD,aes(x,y,lty=branch)) +
    geom_line() +
    theme(legend.position="bottom",
	  #           legend.direction="horizontal",
	  legend.title=element_blank()) +
    scale_linetype_manual(values=c("dashed","solid"))
legs1 <- get_legend(dummyP2)


# extract color legend from plotList and delete legend
legs <- get_legend(pList[[1]])
legs <- plot_grid(legs1,legs,nrow=1)
pls <- lapply(pList,function(x){
	      x <- x %+% guides(color=FALSE) 
	      return(x)})

# add Xlab/Ylab to bottom/left center plots 
pls[["Foxs1"]] <- pls[["Foxs1"]] %+% xlab("Normalized Pseudotime")
pls[["Runx1"]] <- pls[["Runx1"]] %+% ylab("Log-Expression")

#Combine all plots in a single figure
expPlot <- plot_grid(plotlist=pls,ncol=3,labels=c("c","","",
						  "d","","",
						  "e","",""))

expPlot <- plot_grid(expPlot,legs,rel_heights=c(1,0.05),ncol=1)
expPlot <- plot_grid(NULL,expPlot,ncol=1,rel_heights=c(0.3,1))
htmp <- plot_grid(NULL,p0[[4]],ncol=1,rel_heights=c(.275,1),labels=c("a","b"),
		  vjust=c(1.5,0.5))

#Draw branches above heatmap
pb1 <- grid.arrange(pb1)
pb2 <- grid.arrange(pb2)
htmp <- htmp + draw_grob(pb2,-0.02,0.78,.35,.2)
htmp <- htmp + draw_grob(pb1,0.4,0.78,.35,.2)

fullP <- plot_grid(htmp,expPlot,ncol=2,rel_widths=c(1,1))

#close graphics device before plotting
dev.off()
# cairo_pdf("Figure3.pdf",width=16.55,height=13.0575)
fullP
# dev.off()
# 

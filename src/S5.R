# Figure S6
library(plyr)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridGraphics)
library(gridExtra)
library(scran)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# ---- Robustness -----
condComb <- c("NP","G")
set <- c("Bsl-G","Bsl","Myo","Prc")
keepCells <- pD$keep & !(pD$SuperCluster %in% set) & pD$Condition %in% condComb

# Cell filtering
m.vp <- m[,keepCells]
pD.vp <- pD[keepCells,]

# Gene filtering
keep <- rowMeans(m.vp)>0.01
m.vp <- m.vp[keep,]
fD.vp <- fD[keep,]

# Normalize
m.norm <- t(t(m.vp)/pD.vp$sf)

out <- NULL
features <- c("HVG","selected","PCA","all") 
subsampls <- c(1,0.5,0.25)

for (subsampl in subsampls) {
for (feats in features) {
    pD.cur <- pD.vp
    
    if (feats=="HVG") {
    # Highly variable genes 
    var.des <- trendVar(log2(m.norm+1),trend="semiloess")
    var.out <- decomposeVar(log2(m.norm+1),var.des)
    o <- order(var.out$mean)
    hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >=0.5),]

    # Prepare expression matrix
    exps <- m.norm[rownames(hvg.out),]
    exps <- t(log(exps+1))
    }

    if (feats=="selectedGenes") {
    genes <- c("Csn2","Gata3","Prlr","Elf5","Esr1","Pgr","Aldh1a3","Wap",
	       "Tspan8","Krt18","Krt8","Fgfr1","Areg","Fgfr2",
		"Notch1","Notch3","Foxc1","Zeb2")
    exps <- m.norm[fD.vp$symbol %in% genes,]
    exps <- t(log(exps+1))
    }

    if (feats=="PCA") {
    pcs <- prcomp(t(log(m.norm+1)))
    exps  <-pcs$x[,1:50]
    }

    if (feats=="all") {
    exps <- t(log(m.norm+1))
    }

    # Subsampling
    smplsz <- round(subsampl*nrow(exps))
    set.seed(rnd_seed)
    smpl <- sample(rownames(exps),size=smplsz)
    pD.cur <- dplyr::filter(pD.cur, as.character(barcode) %in% smpl)
    exps <- exps[as.character(pD.cur$barcode),]

    # Compute Diffusion map
    set.seed(rnd_seed)
    dm <- DiffusionMap(exps,n_eigs=20,rotate=TRUE)
    pD.cur$DC1 <- eigenvectors(dm)[,1]
    pD.cur$DC2 <- eigenvectors(dm)[,2]
    pD.cur$features <- feats
    pD.cur$SubSample <- paste0(subsampl*100,"%")
    out <- rbind(out,pD.cur)
}
}

# Set color scheme
cols <- levels(pD.vp$Colors)[levels(pD.vp$SubCluster) %in% unique(pD.vp$SubCluster)]

# Plot
out$SubSample <- factor(out$SubSample, levels=c("100%","50%","25%"))
p <- ggplot(out, aes(DC1,DC2, color=SubCluster)) +
    geom_point(size=0.8) +
    scale_color_manual(values=cols) +
    facet_grid(features~SubSample) +
    theme(axis.text=element_blank(),
	  axis.ticks=element_blank(),
	  strip.background=element_blank(),
	  strip.text=element_text(face="bold"),
	  legend.position="bottom",
	  legend.direction="horizontal",
	  ) 

# Run script if not done yet
if (!file.exists("../data/Robjects/secondRun_2500/ExpressionList_Monocle.rds")){
    source("Monocle.R")
}

monoc <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_Monocle.rds")
monoc.plt <- monoc[["plot"]] 
monoc.plt <- monoc[["plot"]] %+% guides(color=FALSE) %+% scale_color_manual(values=cols)

# Combine plots
subp0 <- plot_grid(monoc.plt,NULL,NULL,nrow=1,labels=c("a","",""))
cairo_pdf("../paper/figures/S5.pdf",width=11.69,height=8.27)
plot_grid(subp0,p,ncol=1,labels=c("","b"),rel_heights=c(1,1.5))
dev.off()

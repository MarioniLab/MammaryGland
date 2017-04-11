# Figure S5

library(plyr)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridGraphics)
library(gridExtra)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]


# ---- ScreePlots ----

condComb <- c("NP","G")
sets <- list(NULL,c(6,7,9)) #order important for Robustness section
out <- list()
for (i in seq_along(sets)) {
    set <- sets[[i]]
    keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier & !(pD$cluster %in% set) &
	pD$Condition %in% condComb

    # Cell filtering
    m.vp <- m[,keepCells]
    pD.vp <- pD[keepCells,]

    # Gene filtering
    keep <- rowMeans(m.vp)>0.1
    m.vp <- m.vp[keep,]
    fD.vp <- fD[keep,]

    # Normalize
    m.norm <- t(t(m.vp)/pD.vp$sf)

    # HVG
    brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
    fD.vp$highVar <- fD.vp$id %in% brennecke

    # Scree plot
    exps <- m.norm[fD.vp$highVar,]
    exps <- t(log(exps+1))
    set.seed(rnd_seed)
    dm <- DiffusionMap(exps,n_eigs=20,k=50)
    plot(eigenvalues(dm),pch=20,xlab="Component",ylab="Eigenvalue")
    lines(eigenvalues(dm),pch=20,xlab="Component",ylab="Eigenvalue")
    g <- grab_grob()
    dev.off()
    out[[paste0("g",i)]] <- grid.arrange(g)
}

g1 <- out[["g1"]]
g2 <- out[["g2"]]

# ---- Robustness -----

out <- NULL
features <- c("HVG","selected","PCA","all") 
subsampls <- c(1,0.5,0.25)

for (subsampl in subsampls) {
for (feats in features) {
    pD.cur <- pD.vp
    
    if (feats=="HVG") {
    brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
    exps <- m.norm[brennecke,]
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
    pD.cur <- filter(pD.cur, as.character(barcode) %in% smpl)
    exps <- exps[as.character(pD.cur$barcode),]

    # Compute Diffusion map
    set.seed(rnd_seed)
    dm <- DiffusionMap(exps,n_eigs=20,k=50)
    pD.cur$DC1 <- eigenvectors(dm)[,1]
    pD.cur$DC2 <- eigenvectors(dm)[,2]
    pD.cur$features <- feats
    pD.cur$SubSample <- paste0(subsampl*100,"%")
    out <- rbind(out,pD.cur)
}
}

# Set color scheme
clusts <- as.numeric(as.character(sort(unique(pD.vp$cluster))))
cols <- brewer.pal(n=9,name="Paired")[clusts]
names(cols) <- clusts

# Plot
out$SubSample <- factor(out$SubSample, levels=c("100%","50%","25%"))
p <- ggplot(out, aes(DC1,DC2, color=cluster)) +
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

# Read in Monocle plot
pal <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,8,9)]
monoc <- readRDS("../data/Robjects/ExpressionList_Monocle.rds")
monoc.plt <- monoc[["plot"]] %+% guides(color=FALSE)

# Combine plots
subp0 <- plot_grid(g1,g2,monoc.plt,nrow=1,labels="auto")
cairo_pdf("../paper/figures/S3.pdf",width=11.69,height=8.27)
plot_grid(subp0,p,ncol=1,labels=c("","d"),rel_heights=c(1,1.5))
dev.off()

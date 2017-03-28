library(plyr)
library(splines)
library(destiny)
library(dplyr)
library(reshape)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
source("functions.R")

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

pD$Condition <- mapvalues(pD$Condition, from=c("V","P","L","I"),
			  to=c("Virgin", "Pregnancy",
			       "Lactation", "Post-Involution"))

# Cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier 
m <- m[,keepCells]
pD <- pD[keepCells,]

pD$cluster <- mapvalues(pD$cluster, from=c(1,2,3,4,5,6,7,9,10),
			  to=c(1,2,3,4,5,6,7,8,9))
#Normalize
m.norm <- t(t(m)/pD$sf)
rownames(m.norm) <- fD$symbol


##########################################
#
#
# Diffusion Maps
#
#
##########################################


condComb <- c("Virgin","Pregnancy")

# For Basal and Luminal cells
#
#

keepCells <- pD$Condition %in% condComb
m.vp <- m[,keepCells]
pD.vp <- pD[keepCells,]

#Gene filtering
keep <- rowMeans(m.vp)>0.1
m.vp <- m.vp[keep,]
fD.vp <- fD[keep,]

#
m.norm <- t(t(m.vp)/pD.vp$sf)

brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
fD.vp$highVar <- fD.vp$id %in% brennecke

library(gridGraphics)
library(gridExtra)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

exps <- m.norm[fD.vp$highVar,]
exps <- t(log(exps+1))
set.seed(rnd_seed)
dm <- DiffusionMap(exps,n_eigs=20,k=50)
library(scatterplot3d)
pal <- brewer.pal(n=9,name="Paired")
cols <- mapvalues(pD.vp$cluster,levels(pD.vp$cluster)[-1],
		  pal)
dms <- eigenvectors(dm)[,1:3]
plot(eigenvalues(dm),pch=20,xlab="Component",ylab="Eigenvalue")

g0 <- grab_grob()
g0 <- grid.arrange(g0)

## Extract as grob for plot

#####################################
# Only luminal cells in P and V
#
#####################################

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- !(pD$cluster %in% c(6,7,9)) & pD$Condition %in% condComb
m.vp <- m[,keepCells]
pD.vp <- pD[keepCells,]

#Gene filtering
keep <- rowMeans(m.vp)>0.1
m.vp <- m.vp[keep,]
fD.vp <- fD[keep,]

#
m.norm <- t(t(m.vp)/pD.vp$sf)

brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
fD.vp$highVar <- fD.vp$id %in% brennecke


# Compute Diffusion map
set.seed(rnd_seed)
exps <- m.norm[fD.vp$highVar,]
exps <- t(log(exps+1))
dm <- DiffusionMap(exps,n_eigs=20,k=50)
plot(eigenvalues(dm),pch=20,xlab="Component",ylab="Eigenvalue")
g1 <- grab_grob()
g1 <- grid.arrange(g1)

########### Compute DPT
set.seed(rnd_seed)
dcs <- eigenvectors(dm)[,1:2]
t1 <- which.min(dcs[,2])
t2 <- which.min(dcs[,1])
t3 <- which.max(dcs[,1])
dpt <- DPT(dm, branching=TRUE, tips=c(t1,t2,t3))
root <- which(dpt@tips[,1])[3]
rootdpt <- paste0("DPT",root)
plot(dpt,pch=20,dcs=c(1,2),col_by=rootdpt)

pD.vp$dpt<- dpt[["dpt"]]
########### Compute DPT

out <- NULL
features <- c("HVG","selected","PCA","all") 
subsampls <- c(1,0.5,0.25)

for (subsampl in subsampls) {
for (feats in features) {
    pD.cur <- pD.vp
    if (feats=="HVG") {
    brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
    fD.vp$highVar <- fD.vp$id %in% brennecke
    exps <- m.norm[fD.vp$highVar,]
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
    pD.cur$SubSample <- subsampl
    out <- rbind(out,pD.cur)
}
}

clusts <- as.numeric(as.character(sort(unique(pD.vp$cluster))))
cols <- brewer.pal(n=9,name="Paired")[clusts]
names(cols) <- clusts

library(cowplot)
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


monoc <- readRDS("../data/Robjects/ExpressionList_Monocle.rds")
monoc.plt <- monoc[["plot"]] %+% guides(color=FALSE)
pD.dm <- select(pD.vp, barcode,dpt)
pD.mon <- select(monoc[["pD"]], barcode, Pseudotime)
pD.comp <- inner_join(pD.dm,pD.mon,by="barcode")

corP <- ggplot(pD.comp, aes(x=rank(dpt), y=rank(Pseudotime))) +
    geom_point()

subp0 <- plot_grid(g0,g1,monoc.plt,nrow=1,labels="auto")
# subp1 <- plot_grid(subp0,monoc.plt,nrow=1,rel_widths=c(0.5,1),rel_heights=c(0.5,1),labels=c("","c"))
cairo_pdf("../paper/figures/S5.pdf",width=14.28,height=14.28)
plot_grid(subp0,p,ncol=1,labels=c("","d"),rel_heights=c(1,1.5))
dev.off()

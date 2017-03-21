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

rnd_seed <- 300
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

exps <- m.norm[fD.vp$highVar,]
exps <- t(log(exps+1))
dm <- DiffusionMap(exps,n_eigs=20,k=50)
library(scatterplot3d)
pal <- brewer.pal(n=9,name="Paired")
cols <- mapvalues(pD.vp$cluster,levels(pD.vp$cluster)[-1],
		  pal)
dms <- eigenvectors(dm)[,1:3]

## 3d Plot for Basal cells
scatterplot3d(dms[,2],-dms[,1],dms[,3],color=cols,pch=20,angle=40,
	      cex.symbols=1.5,
	      mar=c(5,3,-0.1,3)+0.1,
	      xlab="DC2",
	      ylab="DC1",
	      zlab="DC3",
	      box=FALSE)

## Extract as grob for plot
library(gridGraphics)
library(gridExtra)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

g <- grab_grob()
g <- grid.arrange(g)

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
# plot(eigenvalues(dm)[1:20])

# Pseudotime and branching
set.seed(rnd_seed)
dcs <- eigenvectors(dm)[,1:2]
t1 <- which.min(dcs[,2])
t2 <- which.min(dcs[,1])
t3 <- which.max(dcs[,1])
dpt <- DPT(dm, branching=TRUE, tips=c(t1,t2,t3))
root <- which(dpt@tips[,1])[3]
rootdpt <- paste0("DPT",root)
# plot(dpt,pch=20,dcs=c(1,2),col_by=rootdpt)

# Rename branches
branch <- dpt@branch[,1]
branch[is.na(branch)]  <- "Decision stage"
branch[branch==1] <- "Root"
branch[branch==2] <- "Secretory lineage"
branch[branch==3] <- "Hormone-sensing lineage"
pD.vp$branch <- factor(branch,levels=c("Root","Decision stage",
				    "Secretory lineage",
				    "Hormone-sensing lineage"))
pD.vp$dpt<- dpt[["dpt"]]
pD.vp$DC1 <- eigenvectors(dm)[,1]
pD.vp$DC2 <- eigenvectors(dm)[,2]

clusts <- as.numeric(as.character(sort(unique(pD.vp$cluster))))
cols <- brewer.pal(n=9,name="Paired")[clusts]
names(cols) <- clusts


## Plots
p0 <- ggplot(pD.vp, aes(x=DC1,y=-DC2, color=cluster)) +
    geom_point(size=2, pch=20) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_color_manual(values=cols)+
    guides(colour=FALSE)

pal <- brewer.pal(n=9,name="Paired")[c(1,2,3,5,6,7,8,9)]
forcLeg <- filter(pD, !cluster %in% c(4))
clustLeg <- ggplot(forcLeg, aes(x=tSNE1,y=tSNE2,color=cluster)) +
    geom_point() +
    scale_color_manual(values=pal) +
    theme(legend.direction="horizontal",
	  legend.title=element_blank(),
	  legend.position="bottom")+
    guides(colour=guide_legend(nrow=1,
			       override.aes=list(size=3)))

clustLeg <- get_legend(clustLeg)
    

p <- ggplot(pD.vp, aes(x=DC1,y=-DC2, color=dpt)) +
    geom_point(size=1, pch=20) +
    scale_color_viridis(option="magma",begin=1,end=0) +
    theme(legend.direction="horizontal")


# Aldehyde trend
ald <- filter(fD.vp, symbol %in% c("Aldh1a3")) %>% .$id
fDC2Trend <- data.frame("DC2"=-pD.vp$DC2,
		   "Aldh1a3"=unname(log2(t(m.norm)[,ald]+1)))
DC2Trend <- ggplot(fDC2Trend, aes(x=DC2, y=Aldh1a3)) +
    geom_point(size=0.5,color="#7570b3",alpha=0.5) +
    geom_smooth(method="loess",color="#7570b3") +
    #     geom_ribbon(aes(ymin = 0,ymax = predict(loess(Aldh1a3 ~ DC2 ))),
    #                     alpha = 0.3,fill = 'dodgerblue')+
#     ggtitle("Progenitor marker expression") +
    ylab("Log-Expression") +
    xlab("") +
    theme(axis.line.y=element_blank(),
	  axis.line.x=element_line(colour="black"),
	  axis.text.y=element_blank(),
	  axis.ticks.y=element_blank()
	  ) +
    coord_flip()

# Csn2 and Esr1 trend
genes <- filter(fD.vp, symbol %in% c("Csn2","Foxa1")) %>% .$id
fDC1Trend <- data.frame("DC1"=pD.vp$DC1,
		   "Csn2"=unname(log2(t(m.norm)[,genes[1]]+1)),
		   "Foxa1"=unname(log2(t(m.norm)[,genes[2]]+1)),
		   "Aldh1a3"=NA)
fDC1Trend <- melt(fDC1Trend,id="DC1",variable_name="Gene")
DC1Trend <- ggplot(fDC1Trend, aes(x=DC1, y=value, color=Gene)) +
    geom_point(size=0.5,alpha=0.5) +
    geom_smooth(method="glm",formula="y~ns(x,df=5)",
		method.args=list(family="quasipoisson")) +
    theme(axis.line.x=element_blank(),
	  axis.line.y=element_line(colour="black"),
	  axis.text.x=element_blank(),
	  axis.ticks.x=element_blank()
	  ) +
    scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3")) +
    ylab("Log-Expression") +
    guides(color=guide_legend(override.aes=list(fill=NA,size=2))) +
    theme(legend.direction="vertical",legend.title=element_blank()) +
    xlab("") 

leg1 <- get_legend(p)
p <- p %+% guides(colour=FALSE)

leg2 <- get_legend(DC1Trend)

DC1Trend <- DC1Trend %+% guides(colour=FALSE)

subP <- plot_grid(p,DC2Trend,DC1Trend,align="hv",rel_widths=c(1,0.6),
		  rel_heights=c(1,0.6))

subP1 <- subP + draw_grob(leg1,0.7,-0.01,1/3,0.5) + draw_grob(leg2,0.7,-0.11,1/3,0.5)
subP1 <- plot_grid(NULL,subP1,NULL,labels=c("","c",""),nrow=1,
		   rel_widths=c(0.2,1,0.2))


subPab <-plot_grid(g,p0,labels=c("a","b",nrow=2)) 
subP0 <- plot_grid(subPab,clustLeg,nrow=2,rel_heights=c(1,0.1))

cairo_pdf("Figure2.pdf",width=12.41,height=15.54)
plot_grid(subP0,subP1,labels=c("",""),nrow=2)
dev.off()

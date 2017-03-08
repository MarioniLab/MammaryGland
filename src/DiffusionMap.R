library(plyr)
library(dplyr)
library(destiny)
library(ggplot2)
source("functions.R")
library(RColorBrewer)
library(rgl)

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier & !(pD$cluster %in% c(6,7,10)) 
m <- m[,keepCells]
pD <- pD[keepCells,]

#Gene filtering
keep <- rowMeans(m)>0.1
m <- m[keep,]
fD <- fD[keep,]

# Genes
#Pre-selected set
# genes <- c("Csn2","Gata3","Prlr","Procr","Elf5","Esr1","Pgr","Aldh1a3","Wap",
#            "Tspan8","Krt18","Krt8","Krt14","Krt5","Fgfr1","Areg","Fgfr2",
#            "Notch1","Notch3","Foxc1","Zeb2")
# fD$highVar <- fD$symbol %in% genes

#Normalize
m.norm <- t(t(m)/pD$sf)

brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)
fD$highVar <- fD$id %in% brennecke

exps <- m.norm[fD$highVar,]
exps <- t(log(exps+1))
dm <- DiffusionMap(exps,n_eigs=20,k=50)
dpt <- DPT(dm, branching=TRUE)
plot(eigenvalues(dm)[1:20])

root <- which(dpt@tips[,1])[3]
rootdpt <- paste0("DPT",root)
plot(dpt,pch=20,dcs=c(1,2),col_by=rootdpt)
plot(dpt,pch=20,dcs=c(1,2),col_by="branch")

# palet <- brewer.pal(length(levels(pD$Condition)),name="Paired")
# cols <- mapvalues(pD$Condition,levels(pD$Condition),palet)
# plot3d(eigenvectors(dm)[,c(1,2,3)],col=cols,radius=0.001, type="s")


#Expression of Genes in the differentiation trajectories
branch <- dpt@branch[,1]
pseudot<- dpt[[rootdpt]]
pseudot[branch==3 & !is.na(branch)] <- -1*pseudot[branch==3 & !is.na(branch)]
gn <- filter(fD, symbol=="Cdc20") %>% .$id
y <- log((t(m.norm)[,gn])+1)
dat <- data.frame(pseudot,y,cluster=pD$cluster,branch=factor(branch),
		  Condition=pD$Condition) %>%
    filter(!is.na(branch) & branch!=3) %>%
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
mat <- m.norm[,branch==1 & !is.na(branch)]
keep <- rowMeans(mat) > 0.5
mat <- mat[keep,]
brennecke <- BrenneckeHVG(mat,suppress.plot=TRUE,fdr=.01)
mat <- log2(mat[brennecke,order(pseudot[branch==1 & !is.na(branch)])]+1)
mat <- mat - rowMeans(mat)
rownames(mat) <- fD$symbol[fD$id %in% brennecke]
pheatmap(mat,
	 cluster_cols=FALSE,
	 #          clustering_distance_rows="correlation",
	 #          clustering_method="ward.D2",
	 show_rownames=TRUE,
	 fontsize=6,
	 show_colnames=FALSE)






pD$branch <- dpt@branch[,1]
pD$DC1 <- eigenvectors(dm)[,1]
pD$DC2 <- eigenvectors(dm)[,2]
p <- ggplot(pD, aes(x=DC1,y=DC2, color=cluster)) +
    geom_point( size=3, pch=20) +
    facet_wrap(~branch) +
    theme_bw() +
    scale_color_brewer(palette="Paired")
p

plot(dpt,dcs=1:3,pch=20,col_by="branch")


#High Var


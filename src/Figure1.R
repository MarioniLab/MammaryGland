#########################################
#
# Figure 1
#
#########################################
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
source("functions.R")

rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]


#Rename Condition for plot
pD$Condition <- mapvalues(pD$Condition, from=c("NP","G","L","PI"),
			  to=c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "14d Post Natural Involution"))

#Remove previously identified outlier and immune cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier 
m <- m[,keepCells]
pD <- pD[keepCells,]

#t-SNE colored by Condition
p0 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=Condition)) +
    geom_point(size=1.5) +
    #     ggtitle("Conditions") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

# t-SNE colored by cluster
p1 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=cluster)) +
    geom_point(size=1.5) +
    scale_color_brewer(palette="Paired")+
    #     ggtitle("Cluster") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

plot_grid(p0,NULL,p1,align="h",nrow=1,rel_widths=c(1,0.2,1))

# ---- Heatmap ----

# Select key genes for each cluster as well as luminal basal markers (general)
general <- c("Acta2","Mylk","Krt5","Krt14","Cnn1","Trp63","Epcam","Krt18","Krt8")
c1 <- c("Cited1","Prlr","Esr1")
c2 <- c("Fabp3","Thrsp","Csn2","Csn1s2a","Glycam1","Olah","Rspo1")
c3 <- c("Foxa1","Aldh1a3","Kit","Cd14")
c4 <- c("Ltf","Hp","Pdk4","C4b","Chil1","Vegfa","Slpi")
c5 <- c("Lypd3","Gpc6")
c6 <- c("1500015O10Rik","Col7a1","Moxd1","Mia","Emid1","Pdpn","Col9a2","Fbln2","Igfbp3","Fst","Il17b","Bmp7")
c7 <- c("Oxtr","Krt15","Lep","Igfbp6")
c8 <- c("Pip","Apod","Prss2","Cbr2","Dusp4")
c9 <- c("Gng11","Procr","Igfbp7","Nrip2","Notch3","Zeb2")

# Combine in order for heatmap
genes <- c(general,c1,c3,c5,c4,c8,c2,c6,c7,c9)

# Subsample cells from large clusters
set.seed(rnd_seed)
subsP <- filter(pD, cluster %in% c(1,2,3,4,5,6)) %>%
    group_by(cluster) %>%
    do(sample_n(.,100))

# Combine with remaining clusters and relevel factor according to order in plot
ord <- filter(pD, cluster %in% c(7,8,9)) %>%
    bind_rows(.,subsP) %>%
    mutate(cluster=factor(cluster,levels=c(1,3,5,4,8,2,6,7,9))) %>%
    arrange(cluster,Condition)

#Normalize data with sizefactors and replace ENSEMBL IDs by gene symbols
m.norm <- t(t(m)/pD$sf)
rownames(m.norm) <- fD$symbol

# Prepare expression matrix for heatmap
mheat <- m.norm[genes,as.character(ord$barcode)]
mheat <- log2(mheat +1)
mheat <- mheat/apply(mheat,1,max) # Scale to 0-1 for visualization

# Prepare Annotation data.frame for heatmap
annoCol <- data.frame("Cluster"=as.factor(ord$cluster),
		      "Stage"=ord$Condition)
rownames(annoCol) <- as.character(ord$barcode)

# Condition color scheme as in F1b
forcol <- ggplot_build(p0)
condColors <- unique(arrange(forcol$data[[1]],group) %>% .[["colour"]])
names(condColors) <- c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "14d Post Natural Involution")
# Cluster color scheme as in F1c
forcol <- ggplot_build(p1)
clustColors <- unique(arrange(forcol$data[[1]],group) %>% .[["colour"]])
clustColors <- clustColors[c(1,3,5,4,8,2,6,7,9)]
names(clustColors) <- c(1,3,5,4,8,2,6,7,9)

# Set Color schemes for annotation data frame 
annoColors <- list("Stage"=condColors,
		   "Cluster"=clustColors)

# Plot heatmap
p <-  pheatmap(mheat,
	 cluster_rows=FALSE,
	 cluster_cols=FALSE,
         show_rownames=TRUE,
         show_colnames=FALSE,
         annotation_legend=FALSE,
	 annotation_col=annoCol,
	 gaps_col=c(100,200,300,400,426,526,626,661),
         gaps_row=9,
	 annotation_colors=annoColors,
	 fontsize=8)


# Combine all plots
subP <- plot_grid(p0,NULL,p1,align="h",nrow=1,rel_widths=c(1,0.2,1),
		  labels=c("b","","c"))
fullP <- plot_grid(subP,p[[4]],nrow=2,
		   labels=c("","d"),vjust=0)
fullPplus0 <- plot_grid(NULL,fullP,nrow=2,
			labels=c("a",""),rel_heights=c(0.3,0.9))

# cairo_pdf("Figure1.pdf",width=12.41,height=17.54)
fullPplus0
# dev.off()

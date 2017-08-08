# Figure 1

library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(gridGraphics)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]


# ---- tSNE-Plots ----

# Rename Condition for plot
pD$Condition <- mapvalues(pD$Condition, from=c("NP","G","L","PI"),
			  to=c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "11d Post Natural Involution"))

# Remove outlier and immune cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

m <- t(t(m)/pD$sf)

# t-SNE colored by Condition
p0 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=Condition)) +
    geom_point(size=1.5) +
    #     ggtitle("Conditions") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

# t-SNE colored by cluster
p1 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=SuperCluster)) +
    geom_point(size=1.5) +
    scale_color_manual(values=levels(pD$SuperColors))+
    #     ggtitle("Cluster") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank())  
#     facet_wrap(~Condition)

plot_grid(p0,NULL,p1,align="h",nrow=1,rel_widths=c(1,0.2,1))

# ---- Dendrogram ----
pD$SubCluster <- factor(pD$SubCluster) #drop unused levels
out <- data.frame(numeric(nrow(m)))
colnames(out) <- levels(pD$SubCluster)[1]
for (clust in levels(pD$SubCluster)) {
    expr <- rowMeans(log2(m[,pD$SubCluster==clust]+1))
    colname <- clust
    out[,colname] <- expr
}

library(dendextend)
dis <- as.dist((1-cor(out,method="spearman"))/2)
dendr <- dis %>% hclust(.,method="ward.D2") %>% as.dendrogram %>%
    set("labels_col",values=c(3,4),k=2) %>%
    #     set("leaves_col",c(3,4)) %>%
    set("leaves_pch",19) %>%
    set("leaves_cex",2) %>%
    set("branches_lwd",3)
par(cex=.9)
plot(dendr,horiz=TRUE,yaxt="n")
p.dendr <- grab_grob()
p.dendr <- grid.arrange(p.dendr)
dev.off()

# ---- Heatmap ----

# Select key genes for each cluster as well as luminal basal markers (general)
general <- c("Acta2","Mylk","Krt5","Krt14","Cnn1","Trp63","Epcam","Krt18","Krt8")
c1 <- c("Cited1","Prlr","Esr1","Areg")
c2 <- c("Rspo1","Atp6v1c2","Fabp3","Thrsp","Wap","Glycam1","Olah")
c3 <- c("Foxa1","Ly6a","Aldh1a3","Kit","Cd14")
# c4 <- c("Ltf","Hp","Pdk4","C4b","Chil1","Vegfa","Slpi")
c5 <- c("Lypd3")
c6 <- c("1500015O10Rik","Col7a1","Moxd1","Mia","Emid1","Pdpn","Col9a2","Fbln2","Igfbp3","Fst","Il17b")
c7 <- c("Oxtr","Krt15","Igfbp6","Igfbp2","Tns1")
# c8 <- c("Pip","Apod")
c9 <- c("Gng11","Procr","Igfbp7","Nrip2","Notch3","Zeb2")

# Combine in order for heatmap
genes <- c(general,c1,c3,c5,c2,c6,c7,c9)

# Subsample cells from large clusters
set.seed(rnd_seed)
ord <- group_by(pD, SuperCluster) %>%
    do(sample_n(.,200)) %>%
    ungroup() %>%
    mutate(SuperCluster=factor(SuperCluster,levels=c("Hsd","Hsp","Lp","Avp","Avd","Bsl","Myo","Prc"))) %>%
    arrange(SuperCluster,Condition) 


#Normalize data with sizefactors and replace ENSEMBL IDs by gene symbols
rownames(m) <- fD$symbol

# Prepare expression matrix for heatmap
mheat <- m[genes,as.character(ord$barcode)]
mheat <- log2(mheat +1)
mheat <- mheat/apply(mheat,1,max) # Scale to 0-1 for visualization

# Prepare Annotation data.frame for heatmap
annoCol <- data.frame("Cluster"=ord$SuperCluster,
		      "Stage"=ord$Condition)
rownames(annoCol) <- as.character(ord$barcode)

# Condition color scheme as in F1b
forcol <- ggplot_build(p0)
condColors <- unique(arrange(forcol$data[[1]],group) %>% .[["colour"]])
names(condColors) <- c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "11d Post Natural Involution")

# Cluster color scheme as in F1c
forcol <- ggplot_build(p1)
clustColors <- unique(arrange(forcol$data[[1]],group) %>% .[["colour"]])
clustColors <- clustColors[as.character(levels(ord$SuperColors))]
clustColors <- as.character(unique(ord$SuperColors))
names(clustColors) <- levels(ord$SuperCluster)

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
	 gaps_col=c(200,400,600,800,1000,1200,1400),
         gaps_row=9,
	 annotation_colors=annoColors,
	 fontsize=8)


# Combine all plots
topRow <- plot_grid(NULL,NULL,p0,labels=c("a","","b"),rel_widths=c(1,0.4,1),nrow=1)
secondRow <- plot_grid(p.dendr,NULL,p1,align="h",nrow=1,rel_widths=c(1,0.4,1),
		  labels=c("c","","d"))
fullP <- plot_grid(topRow,secondRow,p[[4]],nrow=3,
		   labels=c("","","e"),vjust=0)

cairo_pdf("../paper/figures/Figure1.pdf",width=12.41,height=17.54)
fullP
dev.off()

# Figure 2

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
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]


# Remove outlier and immune cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

# Normalize
m <- t(t(m)/pD$sf)


# Rename Condition for plot
pD$Condition <- mapvalues(pD$Condition, from=c("NP","G","L","PI"),
			  to=c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "11d Post Natural Involution"))


# ---- Dendrogram ----

pD$SubClusterNumbers <- factor(pD$SubClusterNumbers) #drop unused levels
out <- data.frame(numeric(nrow(m)))
colnames(out) <- levels(pD$SubClusterNumbers)[1]
for (clust in levels(pD$SubClusterNumbers)) {
    expr <- rowMeans(log2(m[,pD$SubClusterNumbers==clust]+1))
    colname <- clust
    out[,colname] <- expr
}

library(dendextend)
dis <- as.dist((1-cor(out,method="spearman"))/2)
dendr <- dis %>% hclust(.,method="ward.D2") %>% as.dendrogram %>%
    set("labels_col",values=c(3,4),k=2) %>%
    set("leaves_pch",19) %>%
    set("leaves_cex",2) %>%
    set("branches_lwd",3)
par(cex=.6)
plot(dendr,horiz=TRUE,yaxt="n")
p.dendr <- grab_grob()
p.dendr <- grid.arrange(p.dendr)
dev.off()

# ---- Heatmap ----

# Select key genes for each cluster
c1 <- c("Cited1","Prlr","Esr1","Areg")
c2 <- c("Rspo1","Atp6v1c2","Fabp3","Thrsp","Wap","Glycam1","Olah")
c3 <- c("Foxa1","Ly6a","Aldh1a3","Kit","Cd14")
c5 <- c("Lypd3")
c6 <- c("1500015O10Rik","Col7a1","Moxd1","Mia","Emid1","Pdpn","Col9a2","Fbln2","Igfbp3","Fst","Il17b")
c7 <- c("Oxtr","Krt15","Igfbp6","Igfbp2","Tns1")
c9 <- c("Gng11","Procr","Igfbp7","Nrip2","Notch3","Zeb2")

# Combine in order for heatmap
genes <- c(c1,c3,c5,c2,c6,c7,c9)

# Subsample cells from large clusters
set.seed(rnd_seed)
subsP <- filter(pD, !SubCluster %in% c("Hsd-G","Avd-L")) %>%
    group_by(SubCluster) %>%
    do(sample_n(.,100))

# Combine with remaining clusters and relevel factor according to order in plot
ord <- filter(pD, SubCluster %in% c("Hsd-G","Avd-L")) %>%
    bind_rows(.,subsP) %>%
    mutate(SubCluster=factor(SubCluster,levels=c("Hsd-NP","Hsd-PI","Hsd-G","Hsp-NP","Hsp-PI","Lp-NP","Lp-PI","Avp-G",
					   "Avp-L","Avd-G","Avd-L","Bsl","Bsl-G","Myo","Prc"))) %>%
    arrange(SubCluster,Condition)

#Replace ENSEMBL IDs by gene symbols
rownames(m) <- fD$symbol

# Prepare expression matrix for heatmap
mheat <- m[genes,as.character(ord$barcode)]
mheat <- log2(mheat +1)
mheat <- mheat/apply(mheat,1,max) # Scale to 0-1 for visualization

# Prepare Annotation data.frame for heatmap
annoCol <- data.frame("Cluster"=ord$SubCluster,
		      "Stage"=ord$Condition)
rownames(annoCol) <- as.character(ord$barcode)

# ensure to have same color scheme as in figure 1
p0 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=Condition)) +
    geom_point(size=1.5) +
    #     ggtitle("Conditions") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

# Condition color scheme as in F1b
forcol <- ggplot_build(p0)
condColors <- unique(arrange(forcol$data[[1]],group) %>% .[["colour"]])
names(condColors) <- c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "11d Post Natural Involution")

# Cluster color scheme as in F1c
p1 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=SubCluster)) +
    geom_point(size=1.5) +
    scale_color_manual(values=levels(pD$Colors))+
    #     ggtitle("Cluster") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank())  

forcol <- ggplot_build(p1)
clustColors <- unique(arrange(forcol$data[[1]],group) %>% .[["colour"]])
clustColors <- clustColors[as.character(levels(ord$Colors))]
clustColors <- as.character(unique(ord$Colors))
names(clustColors) <- levels(ord$SubCluster)

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
	 gaps_col=c(263,463,663,863,1052,1252,1352),
	 annotation_colors=annoColors,
	 fontsize=9)

fullP <- plot_grid(p.dendr,p[[4]],nrow=2)

# uncomment to save plots
# ggsave(filename="../paper/figures/f2_a.pdf",p[[4]],width=14,height=7)
# ggsave(filename="../paper/figures/f2_b.svg",p.dendr)
# dev.off()
# cairo_pdf("../paper/figures/Figure2.pdf",width=12.41,height=17.54)
# fullP
# dev.off()

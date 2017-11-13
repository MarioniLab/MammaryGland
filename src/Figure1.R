# Figure 1

library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(pheatmap)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# Rename Condition for plot
pD$Condition <- mapvalues(pD$Condition, from=c("NP","G","L","PI"),
			  to=c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "11d Post Natural Involution"))

# Remove outlier and immune cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

# t-SNE colored by Condition
p0 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=Condition)) +
    geom_point(size=1.5) +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

# t-SNE colored by cluster
p1 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=SubClusterNumbers)) +
    geom_point(size=1.5) +
    scale_color_manual(values=levels(pD$Colors))+
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank())  

# Prepare data for expression tSNE
m.norm <- t(t(m)/pD$sf)
rownames(m.norm) <- fD$symbol

titles <- genes <- c("Krt18","Krt5")
expr <- log2(t(m.norm)[,genes]+1)
add <- data.frame(expr,
		  barcode=colnames(m))
colnames(add) <- gsub("X","",colnames(add))
forPlot <- left_join(add,pD[,c("barcode","tSNE1","tSNE2")])
forPlot <- melt(forPlot,id=c("barcode","tSNE1","tSNE2")) %>%
    dplyr::rename(Expression=value)
plots <- list()

# color palette as in heatmap
pal <- colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(200)

# Loop accross genes to create expression plots
for(i in seq_along(genes)) {
    gene <- genes[i]
    fP <- filter(forPlot,variable==gene) %>% arrange(Expression)
    p <- ggplot(fP, aes(x=tSNE1, y=tSNE2, color=Expression)) +
	geom_point(size=1) +
	scale_color_gradientn(colors=pal) +
	ggtitle(titles[i]) +
	theme_void(base_size=14) +
	theme(plot.title=element_text(size=rel(1.05))) 
    plots[[gene]] <- p
}

s0 <- plot_grid(p0,p1,labels=c("b","c"))
s1 <- plot_grid(plotlist=plots)
full <- plot_grid(NULL,s0,s1,nrow=3,labels=c("a","","d"))

# uncomment to save files
# ggsave(filename="../paper/figures/f1_b.pdf",p0)
# ggsave(filename="../paper/figures/f1_c.pdf",p1)
# ggsave(filename="../paper/figures/f1_c1.pdf",plots[[1]])
# ggsave(filename="../paper/figures/f1_c2.pdf",plots[[2]])
# 
# cairo_pdf("../paper/figures/Figure1.pdf", width=8.27, height=11.69)
# full
# dev.off()


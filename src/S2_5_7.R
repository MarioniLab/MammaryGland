# S2 S5 and S7

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


# ---- S2 ----

# Rename Condition for plot
pD$Condition <- mapvalues(pD$Condition, from=c("NP","G","L","PI"),
			  to=c("Nulliparous", "14.5d Gestation",
			       "6d Lactation", "11d Post Natural Involution"))

# Remove previously identified outlier and immune cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

# put in rnd order for plotting
set.seed(rnd_seed)
fp1 <- pD[sample(c(1:nrow(pD)),nrow(pD)),]
#t-SNE colored by SampleID
p1 <- ggplot(fp1, aes(x=tSNE1, y=tSNE2, color=SampleID)) +
    geom_point(size=1) +
    scale_color_brewer(palette="Paired")+
    #     ggtitle("Cluster") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

# Normalize
m.norm <- t(t(m)/pD$sf)
rownames(m.norm) <- fD$symbol

titles <- genes <- c("Krt18","Krt8","Krt5","Krt14","Acta2")
add <- data.frame(log2(t(m.norm)[,genes]+1),
		  barcode=colnames(m))
colnames(add) <- gsub("X","",colnames(add))
forPlot <- left_join(add,pD[,c("barcode","tSNE1","tSNE2")])
forPlot <- melt(forPlot,id=c("barcode","tSNE1","tSNE2")) %>%
    dplyr::rename(Expression=value)
plots <- list()
for(i in seq_along(genes)) {
    gene <- genes[i]
    fP <- filter(forPlot,variable==gene) %>% arrange(Expression)
    p <- ggplot(fP, aes(x=tSNE1, y=tSNE2, color=Expression)) +
	geom_point(size=1) +
	scale_color_viridis(guide_legend(title=gene)) +
	ggtitle(titles[i]) +
	theme_void(base_size=14) +
	theme(plot.title=element_text(size=rel(1.05))) 
    plots[[gene]] <- p
}
cairo_pdf("../paper/figures/S2.pdf",width=11.69,height=8.27)
plot_grid(p1,plotlist=plots,labels=c("a","b","","c","",""))
dev.off()

# ---- S5 ----

# t-SNE colored by SampleID
p2 <- ggplot(fp1, aes(x=tSNE1, y=tSNE2, color=factor(SuperCluster))) +
    geom_point(size=1) +
    scale_color_manual(values=levels(fp1$SuperColor))+
    #     ggtitle("Cluster") +
    theme_void(base_size=12) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

titles <- genes <- c("Aldh1a3","Cd14","Kit")
add <- data.frame(log2(t(m.norm)[,genes]+1),
		  barcode=colnames(m))
colnames(add) <- gsub("X","",colnames(add))
forPlot <- left_join(add,pD[,c("barcode","tSNE1","tSNE2")])
forPlot <- melt(forPlot,id=c("barcode","tSNE1","tSNE2")) %>%
    dplyr::rename(Expression=value)
plots <- list()
for(i in seq_along(genes)) {
    gene <- genes[i]
    fP <- filter(forPlot,variable==gene) %>% arrange(Expression)
    p <- ggplot(fP, aes(x=tSNE1, y=tSNE2, color=Expression)) +
	geom_point(size=1) +
	scale_color_viridis(guide_legend(title=gene)) +
	ggtitle(titles[i]) +
	theme_void(base_size=14) +
	theme(plot.title=element_text(size=rel(1.05))) 
    plots[[gene]] <- p
}
cairo_pdf("../paper/figures/S5.pdf",width=7.79,height=8.27)
plot_grid(p2,plotlist=plots)
dev.off()

# ---- S7 ----
titles <- genes <- c("Csn3","Csn1s2a","Csn1s1","Csn2")
add <- data.frame(log2(t(m.norm)[,genes]+1),
		  barcode=colnames(m))
colnames(add) <- gsub("X","",colnames(add))
forPlot <- left_join(add,pD[,c("barcode","tSNE1","tSNE2")])
forPlot <- melt(forPlot,id=c("barcode","tSNE1","tSNE2")) %>%
    dplyr::rename(Expression=value)
plots <- list()
for(i in seq_along(genes)) {
    gene <- genes[i]
    fP <- filter(forPlot,variable==gene) %>% arrange(Expression)
    p <- ggplot(fP, aes(x=tSNE1, y=tSNE2, color=Expression)) +
	geom_point(size=1) +
	scale_color_viridis(guide_legend(title=gene)) +
	ggtitle(titles[i]) +
	theme_void(base_size=14) +
	theme(plot.title=element_text(size=rel(1.05))) 
    plots[[gene]] <- p
}
cairo_pdf("../paper/figures/S7.pdf",width=11.69,height=8.27)
plot_grid(p2,plotlist=plots)
dev.off()

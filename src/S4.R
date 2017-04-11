# S4

source("Figure3.R")
genes1 <- res1$Gene
genes2 <- res2$Gene

#Rename the alveolar cells, so that there are no duplicate names for rows
colnames(m.both) <- c(paste0("Alv.",colnames(m.both)[1:ncol(m.alv)]),
		       colnames(m.both)[(ncol(m.alv)+1):ncol(m.both)])

m.heat1 <- m.both[genes1,]
m.heat2 <- m.both[genes2,]

##Plot
p0 <- pheatmap(m.heat1,
	 cluster_cols=FALSE,
	 cluster_rows=TRUE,
	 clustering_distance_rows="correlation",
	 annotation=annoCol,
	 clustering_method="average",
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 treeheight_row=0,
	 legend=FALSE,
	 annotation_legend=FALSE,
	 gaps_col=ncol(m.alv),
	 show_rownames=FALSE)

p1 <- pheatmap(m.heat2,
	 cluster_cols=FALSE,
	 cluster_rows=TRUE,
	 clustering_distance_rows="correlation",
	 annotation=annoCol,
	 clustering_method="average",
	 show_colnames=FALSE,
	 annotation_colors=annoColors,
	 treeheight_row=0,
	 annotation_legend=TRUE,
	 gaps_col=ncol(m.alv),
	 show_rownames=FALSE)
png("../paper/figures/S4.png",height=1248,width=1248)
comb <- plot_grid(p0[[4]],NULL,p1[[4]],nrow=1,vjust=0.5,rel_widths=c(1,.1,1))
comb
dev.off()

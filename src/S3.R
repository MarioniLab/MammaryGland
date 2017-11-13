# S4 and S7
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(reshape2)
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


# ---- S4 ----

pD$SampleID <- factor(pD$SampleID) #drop unused levels
out <- data.frame(numeric(nrow(m)))
colnames(out) <- levels(pD$SampleID)[1]
for (sample in levels(pD$SampleID)) {
    expr <- rowMeans(log2(m[,pD$SampleID==sample]+1))
    colname <- sample
    out[,colname] <- expr
}
plot(out[,1]~out[,2],pch=20,xlab="NP1",ylab="NP2")
abline(0,1,col="dodgerblue",lwd=2)
legend(x="topleft",legend=paste("Cor=",round(cor(out[,1],out[,2]),digits=4)))
pNP <- grab_grob()
pNP <- grid.arrange(pNP)
dev.off()

plot(out[,3]~out[,4],pch=20,xlab="G1",ylab="G2")
abline(0,1,col="dodgerblue",lwd=2)
legend(x="topleft",legend=paste("Cor=",round(cor(out[,3],out[,4]),digits=4)))
pG <- grab_grob()
pG <- grid.arrange(pG)
dev.off()

plot(out[,5]~out[,6],pch=20,xlab="L1",ylab="L2")
abline(0,1,col="dodgerblue",lwd=2)
legend(x="topleft",legend=paste("Cor=",round(cor(out[,5],out[,6]),digits=4)))
pL <- grab_grob()
pL <- grid.arrange(pL)
dev.off()

plot(out[,7]~out[,8],pch=20,xlab="PI1",ylab="PI2")
abline(0,1,col="dodgerblue",lwd=2)
legend(x="topleft",legend=paste("Cor=",round(cor(out[,7],out[,8]),digits=4)))
pPI <- grab_grob()
pPI <- grid.arrange(pPI)
dev.off()

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

corPlots <- plot_grid(pNP,pG,pL,pPI)
cairo_pdf("../paper/figures/S3.pdf",width=11.69,height=8.27)
plot_grid(p1,corPlots)
dev.off()

# Figure 4

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(reshape2)
library(RColorBrewer)
source("functions.R")
# 
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# ---- LpVolcanoPlot ----

topTab <- read.csv(file="../data/Robjects/Lp_NPvsPI.csv")

# Highlight genes with lactation/immune annotation
lac <- read.csv("../data/Robjects/LactationGenes.csv", stringsAsFactors=FALSE)$x
immuno <- read.csv("../data/Robjects/ImmuneGenes.csv", stringsAsFactors=FALSE)$x

# Add some of the Csn genes that are not annotated in GO
lac <- c(lac,"Csn1s1","Csn1s2a","Lalba","Btn1a1")

# Highlight top DE genes in Volcano plot
topUp <- filter(topTab, FDR < 0.01 & logFC > 0) %>%
    arrange(desc(logFC)) %>% .$symbol %>% as.character() %>% .[1:5]
topDown <- filter(topTab, FDR < 0.01 & logFC < 0) %>%
    arrange(logFC) %>% .$symbol %>% as.character() %>% .[1:5]

# Genes to highlight
interest <- filter(topTab, symbol %in% c(topUp,topDown)) 
imminterest <- filter(topTab, symbol %in% immuno)
lacinterest <- filter(topTab, symbol %in% lac)

# VolcanoPlot
volcano <- ggplot(topTab,aes(x=logFC,y=-log10(FDR))) +
    geom_point(size=2,color="grey50",pch=20) +
    geom_hline(yintercept=2,lty="dashed") +
    geom_vline(xintercept=1,lty="dashed") +
    geom_vline(xintercept=-1,lty="dashed") +
    geom_point(data=interest, aes(x=logFC, y=-log10(FDR)), size=3, color="black", pch=20) +
    geom_point(data=lacinterest, aes(x=logFC, y=-log10(FDR)), size=3, color="dodgerblue", pch=20) +
    geom_point(data=imminterest, aes(x=logFC, y=-log10(FDR)), size=3, color="coral", pch=20) +
    geom_text_repel(data=interest, aes(x=logFC, y=-log10(FDR), label=symbol, fontface="italic")) +
    xlab("Log2 Fold Change") +
    ylab("-Log10(P value)") 

# ---- ProgenitorvsLuminal ----

progenitorDE <- read.csv("../data/Robjects/ProgenitorDE.csv")
# genes to highlight
interest <- filter(progenitorDE, Gene %in% c("Aldh1a3", "Lypd3", "Prlr", "Esr1", "Pgr"))

# fold-change plot
p2 <- ggplot(progenitorDE, aes(x=NullParFC,y=ParousFC)) +
    geom_point(color="grey50",size=2) +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), color="black") +
    geom_text_repel(data=interest, aes(x=NullParFC,y=ParousFC,label=Gene,fontface="italic")) +
    xlab("Lp-NP vs. luminal cells [log2 FC]") +
    ylab("Lp-PI vs. luminal cells [log2 FC]") +
    geom_hline(yintercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    coord_equal(xlim=c(-5,5),ylim=c(-5,5))

# ---- DEstrongerInProgenitor ----

# 
m <- m[fD$keep,pD$keep]
pD <- pD[pD$keep,]
fD <- fD[fD$keep,]

# 
m <- t(t(m)/pD$sf)

rownames(m) <- as.character(fD$symbol)
genes1 <- c("Csn1s1","Lalba",
	   "Lipa",
	   "Cidea",
	   "Xdh","Cd36")
genes2 <- c("Ctsc","Cd14","Tgfb3","Mfge8","Hp","Spp1")
geneList <- list(genes1,genes2)
out <- list()

# loop over the two gene list to create expression plots
for (i in c(1,2)) {
    genes <- geneList[[i]]

    #Create DF
    forPl <- data.frame(t(m)[,c(genes)]+1,
	    barcode=colnames(m))
    add <- select(pD, barcode, SubCluster, SuperCluster, Condition, Colors) %>%
	mutate(barcode=as.character(barcode))
    forPl <- left_join(forPl,add,by="barcode") %>%
	filter(SubCluster %in% c("Hsd-NP","Hsd-PI","Hsp-NP","Hsp-PI","Lp-NP","Lp-PI"))

    forPl <- melt(forPl,id=c("barcode","SuperCluster","SubCluster","Condition","Colors"))

    #color scheme
    cols <- levels(forPl$Colors)[levels(forPl$SubCluster) %in% unique(forPl$SubCluster)]
    names(cols) <- as.character(levels(forPl$SubCluster))[levels(forPl$SubCluster) %in% unique(forPl$SubCluster)]
    forPl$SubCluster <- factor(forPl$SubCluster,levels=c("Lp-PI","Lp-NP","Hsp-PI","Hsp-NP","Hsd-PI","Hsd-NP"))
    cols <- cols[as.character(levels(forPl$SubCluster))]

    #Plot
    out[[i]] <- ggplot(forPl, aes(y=value,x=SubCluster,color=SubCluster)) +
	geom_jitter(size=0.9) +
	stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
		     geom = "crossbar", width = 1,color="black") +
	facet_wrap(~variable,scales="free_y",nrow=1) +
	scale_color_manual(values=cols) +
	ylab("") +
	xlab("") +
	theme(strip.background=element_blank(),
	      strip.text=element_text(face="bold"),
	      legend.position="bottom",
	      legend.direction="horizontal",
	      legend.title=element_blank(),
	      axis.text.y = element_text(size=4),
      	      axis.text.x = element_blank(), 
	      axis.ticks.x = element_blank(),
	      strip.text.x=element_text(face="italic"))+
	guides(colour = FALSE) +
	scale_y_log10()
}
plot(out[[1]])
plot(out[[2]])
dev.off()

#dummy legend
dum <- data.frame(x=c(1,1),y=c(1,1),grp=c("Immune response","Lactation"))
dumleg <- ggplot(dum, aes(x,y,color=grp)) +
    geom_point() +
    scale_color_manual(values=c("coral","dodgerblue")) +
    theme(legend.position="bottom",
	  legend.title=element_blank(),
	  legend.direction="horizontal") +
    guides(colour = guide_legend(override.aes = list(size=3))) 
dumleg <- get_legend(dumleg)


# Combine plots
leg <- plot_grid(NULL,dumleg,nrow=1)
subp <- plot_grid(p2,volcano,nrow=1,labels=c("a","b"),align="h")
subp <- plot_grid(subp,leg,rel_heights=c(1,0.1),nrow=2)

# ---- Lp post-involution is biased towards the alveolar fate ----
dms1 <- read.csv("../data/Robjects/dm_luminal.csv")
dms2 <- read.csv("../data/Robjects/dm_luminal_PI.csv")
dms <- rbind(dms1[,-c(6,7)],dms2)

pD <- right_join(pD,dms,by="barcode")
levels(pD$SubCluster) <- c(levels(pD$SubCluster),"NP and G")
pD$SubCluster[pD$Condition %in% c("NP","G")] <- "NP and G"
pD$SubCluster <- factor(as.character(pD$SubCluster), levels=unique(as.character(pD$SubCluster)))
cols <- unique(as.character(pD$Colors))[c(5,7,3,1,6,8,4,9,2)]
cols[9] <- "#D3D3D3"

# Plot with post-involution cells highlighted
p3 <- ggplot(filter(pD, Condition %in% c("L","PI")), aes(DC1,DC2,color=SubCluster)) +
    geom_point(size=2) +
    xlab("Component 1") +
    ylab("Component 2") +
    geom_point(data=pD,aes(DC1,DC2,color=SubCluster),size=0.5) +
    scale_color_manual(values=cols) +
    theme(legend.position="bottom",
	  legend.direction="horizontal",
	  legend.title=element_blank())

legs <- get_legend(p3)
p3 <- p3 %+% guides(color=FALSE)
ExpPlot  <- plot_grid(plotlist=out,nrow=2,labels=c("c","d"))
subp1 <- plot_grid(p3,legs,nrow=1,labels=c("e",""))
fullP <- plot_grid(subp,ExpPlot,subp1,nrow=3)

# uncomment to save
# cairo_pdf("../paper/figures/Figure5.pdf",width=10.75,height=15.19)
# plot_grid(fullP,NULL,nrow=2,rel_heights=c(1,0.25))
# dev.off()

# Figure 4

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(reshape2)
library(RColorBrewer)
source("functions.R")
# 
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# ---- ImpactStrongestOnProgenitor ----

out <- data.frame()
for (clust in c("C1","C3","C5")) {
    filen <- sprintf("../data/Robjects/secondRun_2500/%s_NPvsPI.csv",clust)
    tmp <- read.csv(file=filen)
    tmp$Cluster <- clust
    out <- rbind(out,tmp)
}
out <- group_by(out, Cluster) %>%
    summarize(Up=sum(FDR < 0.01 & logFC > 0),
	      Down=sum(FDR < 0.01 & logFC < 0)) %>%
    melt()

out$Cluster <- factor(out$Cluster,levels=c("C5","C3","C1"))
p1 <- ggplot(out, aes(x=Cluster,y=value,fill=variable)) +
    geom_bar(stat="identity",position="dodge") +
    coord_flip() +
    ylab("Number of Genes DE") +
    labs(fill="")

# ---- C5DE ----

topTab <- read.csv(file="../data/Robjects/secondRun_2500/C5_NPvsPI.csv")

# Highlight genes with lactation/immune annotation
lac <- read.csv("../data/Robjects/secondRun_2500/LactationGenes.csv", stringsAsFactors=FALSE)$x
immuno <- read.csv("../data/Robjects/secondRun_2500/ImmuneGenes.csv", stringsAsFactors=FALSE)$x

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
    geom_label_repel(data=interest, aes(x=logFC, y=-log10(FDR), label=symbol)) +
    xlab("Log2 Fold Change") +
    ylab("-Log10(P value)") 

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

for (i in c(1,2)) {
    genes <- geneList[[i]]
    #Create DF
    forPl <- data.frame(t(m)[,c(genes)]+1,
	    barcode=colnames(m))
    add <- select(pD, barcode, SubCluster, SuperCluster, Condition) %>%
	mutate(barcode=as.character(barcode))
    forPl <- left_join(forPl,add,by="barcode") %>%
	filter(SubCluster %in% c("C1-NP","C1-PI","C3-NP","C3-PI","C5-NP","C5-PI"))

    forPl <- melt(forPl,id=c("barcode","SuperCluster","SubCluster","Condition"))

    #Plot
    forPl$SubCluster <- factor(forPl$SubCluster,levels=c("C5-NP","C3-NP","C1-NP","C5-PI","C3-PI","C1-PI"))
    out[[i]] <- ggplot(forPl, aes(y=value,x=SubCluster,color=SuperCluster)) +
	geom_jitter(size=0.9) +
	stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
		     geom = "crossbar", width = 1,color="black") +
	facet_grid(~variable) +
	ylab("Expression") +
	xlab("") +
	theme(strip.background=element_blank(),
	      strip.text=element_text(face="bold"),
	      legend.position="bottom",
	      legend.direction="horizontal",
	      axis.text.x = element_blank(),
      	      axis.text.y = element_blank(),
	      axis.ticks = element_blank()) +
	guides(colour = guide_legend(override.aes = list(size=3))) +
	scale_y_log10(breaks=c(2,5,10,25,50,100,200,300,500,1000))
}


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
leg <- plot_grid(dumleg,get_legend(out[[1]]),nrow=1)
out[[1]] <- out[[1]] %+% guides(color=FALSE)
out[[2]] <- out[[2]] %+% guides(color=FALSE)
ExpPlot  <- plot_grid(plotlist=out,nrow=2)
subp <- plot_grid(volcano,ExpPlot,nrow=1,labels=c("a","b"))
subp <- plot_grid(subp,leg,rel_heights=c(1,0.1),nrow=2)

# ---- ProgenitorvsLuminal ----
progenitorDE <- read.csv("../data/Robjects/secondRun_2500/ProgenitorDE.csv")
# genes to highlight
interest <- filter(progenitorDE, Gene %in% c("Aldh1a3","Lypd3",
					"Prlr","Esr1","Pgr"))

# FC plot
p2 <- ggplot(progenitorDE, aes(x=NullParFC,y=ParousFC)) +
    geom_point(color="grey50",size=2) +
    geom_point(data=interest, aes(x=NullParFC, y=ParousFC), color="black") +
    geom_label_repel(data=interest, aes(x=NullParFC,y=ParousFC,label=Gene)) +
    xlab("LFC of C5 vs. luminal cells") +
    ylab("LFC of C4 vs. luminal cells") +
    geom_hline(yintercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    coord_equal(xlim=c(-5,5),ylim=c(-5,5))

# ---- C5 post-involution is biased towards the alveolar fate ----
dms1 <- read.csv("../data/Robjects/secondRun_2500/dm_luminal.csv")
dms2 <- read.csv("../data/Robjects/secondRun_2500/dm_luminal_PI.csv")
dms <- rbind(dms1[,-c(6,7)],dms2)

pD <- right_join(pD,dms,by="barcode")
pD$SuperCluster[pD$Condition %in% c("NP","G")] <- NA


p3 <- ggplot(pD, aes(DC1,DC2,color=SuperCluster)) +
    geom_point(size=2) +
    xlab("Component 1") +
    ylab("Component 2") +
    theme(legend.position="bottom",
	  legend.direction="horizontal")
    #     facet_wrap(~Condition)

subp1 <- plot_grid(p2,p3,nrow=1,labels=c("c","d"))
fullP <- plot_grid(subp,subp1,nrow=2)

cairo_pdf("../paper/figures/Figure4.pdf",width=10.75,height=15.19)
plot_grid(fullP,NULL,nrow=2,rel_heights=c(1,0.75))
dev.off()

# GO-Analysis
library(ggplot2)
library(cowplot)
library("topGO")
library(org.Mm.eg.db)
library(dplyr)
library(edgeR)
library(reshape2)
source("functions.R")

# ---- ImpactStrongestOnProgenitor ----

out <- data.frame()
for (clust in c("Hsd","Hsp","Lp")) {
    filen <- sprintf("../data/Robjects/secondRun_2500/%s_NPvsPI.csv",clust)
    tmp <- read.csv(file=filen)
    tmp$Cluster <- clust
    out <- rbind(out,tmp)
}
out <- group_by(out, Cluster) %>%
    summarize(Up=sum(FDR < 0.01 & logFC > 0),
	      Down=sum(FDR < 0.01 & logFC < 0)) %>%
    melt()

out$Cluster <- factor(out$Cluster,levels=c("Lp","Hsp","Hsd"))
p1 <- ggplot(out, aes(x=Cluster,y=value,fill=variable)) +
    geom_bar(stat="identity",position="dodge") +
    coord_flip() +
    ylab("Number of Genes DE") +
    labs(fill="")


# ---- DEAnalysis4vs5 ----
output <- data.frame()
for (clust in c("Hsd","Hsp","Lp")) { # Lp has to be last for rest of script
    filen <- sprintf("../data/Robjects/secondRun_2500/%s_NPvsPI.csv",clust)
    topTab <- read.csv(file=filen)

    deGenes <- topTab[topTab$FDR < 0.01 & topTab$logFC > 0,"symbol"] 
    allGenes <- topTab$symbol

    # ---- Data ----
    univrs <- allGenes
    alG <- factor(as.numeric(univrs %in% deGenes))
    names(alG) <- univrs

    # ---- GOanalysis ----

    # prepare Data for topGO
    GO.data <- new("topGOdata", description="Lib GO",ontology="BP", allGenes=alG, 
		   annot=annFUN.org, mapping="org.Mm.eg.db",
		   nodeSize=20, ID="symbol")
    result.classic <- runTest(GO.data, statistic="fisher")
    tmp <- GenTable(GO.data, Fisher.classic=result.classic, orderBy="topgoFisher", topNodes=50, numChar=30)

    tmp$Term <- factor(tmp$Term, levels=unique(rev(tmp$Term)))
    tmp$Cluster <- clust
    output <- rbind(tmp,output)
}

output <- output[output$Fisher.classic < 0.01,]

ordr <- group_by(output, Cluster) %>%
    arrange(-log10(as.numeric(Fisher.classic))) %>%
    .$Term

output$Term <- factor(output$Term,levels=unique(ordr))

## Extract gene lists for volcano plot
# better hardcode go-terms at one point to define gene lists
sub.output <- output[output$Cluster=="Lp",]
imterm <- sub.output[as.character(c(5,11,12,21,23,29,37,41,42,44,49,50)),1]
lterm <- sub.output[as.character(c(35,25,6,20)),1]
im.genes <- intersect(unlist(genesInTerm(GO.data,imterm)),deGenes)
la.genes <- intersect(unlist(genesInTerm(GO.data,lterm)),deGenes)
write.csv(im.genes, "../data/Robjects/secondRun_2500/ImmuneGenes.csv", row.names=FALSE, quote=FALSE)
write.csv(la.genes, "../data/Robjects/secondRun_2500/LactationGenes.csv", row.names=FALSE, quote=FALSE)

p2 <- ggplot(output, aes(x=Term, y=-log10(as.numeric(Fisher.classic)))) +
    geom_bar(stat="identity",color="black",fill="white") +
    coord_flip() +
    ylab("-log10(P)") +
    geom_hline(yintercept=2,lty="dashed") +
    xlab("GO-Term [BP]") +
    theme(axis.text.y=element_text(size=4)) +
    facet_grid(Cluster~.,scales="free")

# ---- PCAPlot ----
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# remove cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

# DE Test on NP cells progenitor versus rest
pD.sub <- filter(pD,Condition %in% "NP", !(SubCluster %in% c("Bsl","Bsl-G1","Myo","Prc")))
m.sub <- m[,as.character(pD.sub$barcode)]
keep <- rowMeans(m.sub) > 0.01
m.sub <- m.sub[keep,]
fD.sub <- fD[keep,]
rownames(m.sub) <- fD.sub$symbol

nf <- log(pD.sub$sf/pD.sub$UmiSums)
pD.sub$nf <- exp(nf-mean(nf))
y <- DGEList(counts=m.sub,
	     samples=pD.sub,
	     genes=fD.sub,
	     norm.factors=pD.sub$nf)

choice <- "Lp-NP"
cluster <- factor(as.numeric(pD.sub$SubCluster==choice))
de.design <- model.matrix(~cluster)
y <- estimateDisp(y, de.design, prior.df=0,trend="none")
fit <- glmFit(y, de.design)
res <- glmTreat(fit,lfc=1)
resTab <- topTags(res,n=Inf,sort.by="PValue")
tabNulPar <- resTab$table[1:500,]

# ---- PCA ----
pD.sub <- filter(pD,Condition %in% "PI", !(SubCluster %in% c("Bsl","Bsl-G1","Myo","Prc")))
m.sub <- m[,as.character(pD.sub$barcode)]

genes <- tabNulPar$id %in% rownames(m.sub)
m.sub <- t(t(m.sub)/pD.sub$sf)
m.sub <- m.sub[genes,]
m.sub <- t(log2(m.sub+1))
pc <- prcomp(m.sub)
pD.sub$PC1 <- pc$x[,1]
pD.sub$PC2 <- pc$x[,2]

pcplot <- ggplot(pD.sub, aes(PC1,PC2,color=SuperCluster)) +
    geom_point() 
    #     scale_color_manual(values=pal)

loadngs <- pc$rotation[,1]
loddf <- data.frame("id"=names(loadngs),
		    "Loadings"=loadngs)
ld <- plyr::join(tabNulPar,loddf) %>%
    mutate("Genes"=ifelse(logFC >0,"High in Progenitor","High in Differentiated"))
pLoad <- ggplot(ld, aes(x=Genes,y=Loadings)) +
    geom_boxplot() +
    ylab("PC1 Loadings")

################
subp0 <- plot_grid(p1,p2,labels="auto")
subp1 <- plot_grid(pcplot, pLoad, rel_widths=c(1,.8),labels=c("c","d")) 

fullP <- plot_grid(subp0,subp1,nrow=2)
cairo_pdf("../paper/figures/S6.pdf",width=8.27,height=11.69)
fullP
dev.off()

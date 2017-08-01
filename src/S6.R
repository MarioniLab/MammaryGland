# GO-Analysis
library(ggplot2)
library(cowplot)
library("topGO")
library(org.Mm.eg.db)
source("functions.R")

# ---- DEAnalysis4vs5 ----
output <- data.frame()
for (clust in c("C1","C3","C5")) { # C5 has to be last for rest of script
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
		   nodeSize=10, ID="symbol")
    result.classic <- runTest(GO.data, statistic="fisher")
    tmp <- GenTable(GO.data, Fisher.classic=result.classic, orderBy="topgoFisher", topNodes=50, numChar=10000)

    tmp$Term <- factor(tmp$Term, levels=rev(tmp$Term))
    tmp$Cluster <- clust
    output <- rbind(tmp,output)
}

output <- output[output$Fisher.classic < 0.01,]

library(dplyr)

ordr <- group_by(output, Cluster) %>%
    arrange(-log10(as.numeric(Fisher.classic))) %>%
    .$Term

output$Term <- factor(output$Term,levels=unique(ordr))

p1 <- ggplot(output, aes(x=Term, y=-log10(as.numeric(Fisher.classic)))) +
    geom_bar(stat="identity",color="black",fill="white") +
    coord_flip() +
    ylab("-log10(P)") +
    geom_hline(yintercept=2,lty="dashed") +
    xlab("GO-Term [BP]") +
    facet_wrap(~Cluster,scales="free")



subp1 <- plot_grid(subp0,p1,nrow=2,labels=c("","c"))
# cairo_pdf("../paper/figures/S6.pdf",width=9.92,height=14.028)
plot_grid(subp1,nrow=2,labels=c("","d"))
# dev.off()


## Extract gene lists for volcano plot
# better hardcode go-terms at one point to define gene lists
sub.output <- output[output$Cluster=="C5",]
imterm <- sub.output[as.character(c(5,11,12,21,23,29,37,41,42,44,49,50)),1]
lterm <- sub.output[as.character(c(35,25,6,20)),1]
im.genes <- intersect(unlist(genesInTerm(GO.data,imterm)),deGenes)
la.genes <- intersect(unlist(genesInTerm(GO.data,lterm)),deGenes)
write.csv(im.genes, "../data/Robjects/secondRun_2500/ImmuneGenes.csv", row.names=FALSE, quote=FALSE)
write.csv(la.genes, "../data/Robjects/secondRun_2500/LactationGenes.csv", row.names=FALSE, quote=FALSE)


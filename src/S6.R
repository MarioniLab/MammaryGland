# GO-Analysis
library(ggplot2)
library(cowplot)
library("topGO")
library(org.Mm.eg.db)

# ---- DEAnalysis4vs5 ----
output <- data.frame()
input <- read.csv("../data/Robjects/secondRun_2500/BranchDECluster.R",
		  stringsAsFactor=FALSE)
input$BranchCluster <- paste(input$Branch,input$Cluster,sep="-")

for (clust in unique(input$BranchCluster)) { 

    # DeGenes
    deGenes <- input[input$BranchCluster==clust,"Gene"]

    # Gene universe
    if (substr(clust,1,3)=="Hrm") i <- 1 else i <- 2
    m.smooth <- readRDS("../data/Robjects/secondRun_2500/BranchDEList.rds")[[i]][["mSmooth"]]
    univrs <- rownames(m.smooth[rowMeans(m.smooth)>0.1,])

    # ---- Data ----
    alG <- factor(as.numeric(univrs %in% deGenes))
    names(alG) <- univrs

    # ---- GOanalysis ----

    # prepare Data for topGO
    GO.data <- new("topGOdata", description="Lib GO",ontology="BP", allGenes=alG, 
		   annot=annFUN.org, mapping="org.Mm.eg.db",
		   nodeSize=20, ID="symbol")
    result.classic <- runTest(GO.data, statistic="fisher")
    tmp <- GenTable(GO.data, Fisher.classic=result.classic, orderBy="topgoFisher", topNodes=50, numChar=300)

    tmp$Term <- factor(tmp$Term, levels=unique(rev(tmp$Term)))
    tmp$Cluster <- clust
    output <- rbind(tmp,output)
}

output <- output[output$Fisher.classic < 0.01,]

library(dplyr)
# ordr <- group_by(output, Cluster) %>%
#     arrange(-log10(as.numeric(Fisher.classic))) %>%
#     .$Term

# output$Term <- factor(output$Term,levels=unique(ordr))
output1 <- output[grepl("Alv",output$Cluster),]
output2 <- output[grepl("Hrm",output$Cluster),]

p1 <- ggplot(output1, aes(x=Term, y=-log10(as.numeric(Fisher.classic)))) +
    geom_bar(stat="identity",color="black",fill="white") +
    coord_flip() +
    ylab("-log10(P)") +
    geom_hline(yintercept=2,lty="dashed") +
    xlab("GO-Term [BP]") +
    theme(axis.text.y=element_text(size=7)) +
    facet_grid(Cluster~.,scales="free")

p2 <- ggplot(output2, aes(x=Term, y=-log10(as.numeric(Fisher.classic)))) +
    geom_bar(stat="identity",color="black",fill="white") +
    coord_flip() +
    ylab("-log10(P)") +
    geom_hline(yintercept=2,lty="dashed") +
    xlab("GO-Term [BP]") +
    theme(axis.text.y=element_text(size=7)) +
    facet_grid(Cluster~.,scales="free")

cairo_pdf("../paper/figures/S6.pdf",width=18,height=10)
plot_grid(p1,p2,labels="auto")
dev.off()

# Figure 4

library(dplyr)
library(edgeR)
source("functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

# Loop over all luminal clusters
comps <- list(Hsd=c("Hsd-PI","Hsd-NP"),
	      Hsp=c("Hsp-PI","Hsp-NP"),
	      Lp=c("Lp-PI","Lp-NP"),
	      Bsl=c("Bsl-G","Bsl"))

for (cname in c("Hsd","Hsp","Lp","Bsl")) {
    comp <- comps[[cname]]
    pD.sub <- pD[pD$SubCluster %in% comp,]
    m.sub <- m[,as.character(pD.sub$barcode)]
    keep <- rowMeans(m.sub) > 0.1
    m.sub <- m.sub[keep,]
    fD.sub <- fD[keep,]
    rownames(m.sub) <- fD.sub$symbol

    # DGEList object
    nf <- log(pD.sub$sf/pD.sub$UmiSums)
    pD.sub$nf <- exp(nf-mean(nf))
    y <- DGEList(counts=m.sub,
		 samples=pD.sub,
		 genes=fD.sub,
		 norm.factors=pD.sub$nf)

    # DE
    choice <- comp[1] 
    cluster <- factor(as.numeric(pD.sub$SubCluster==choice))
    de.design <- model.matrix(~cluster)
    y <- estimateDisp(y, de.design, prior.df=0, trend="none")
    fit <- glmFit(y, de.design)
    res <- glmTreat(fit, lfc=1)
    resTab <- topTags(res, n=Inf, sort.by="PValue")
    topTab <- resTab$table

    # Write DE table for supps
    topTab <- select(topTab, id, symbol, logFC, unshrunk.logFC, logCPM, PValue, FDR)
    fname <- sprintf("../data/Robjects/%s_NPvsPI.csv",cname)
    if (cname=="Bsl") {fname <- sprintf("../data/Robjects/%s_NPvsG.csv",cname)}
    write.csv(topTab,file=fname, row.names=FALSE)
}

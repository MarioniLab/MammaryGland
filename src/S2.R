library(scran)
library(plyr)
library(dplyr)
library(knitr)
library(ggplot2)
library(cowplot)
library(Rtsne)
source("functions.R")
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
#relevel pD Condition to get it in logical sequence
rnd_seed <- 300

pD$Condition <- mapvalues(pD$Condition, from=c("V","P","L","I"),
			  to=c("NP", "G",
			       "L", "PI"))
pD$Condition <- factor(pD$Condition, levels=c("NP","G","L","PI"))
pD$SampleID <- mapvalues(pD$SampleID, from=c("V1","V2","P1","P2","L1","L2","I1","I2"),
			  to=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2"))
pD$SampleID <- factor(pD$SampleID,levels=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2")) 


# A more or less fair representation of the lactation samples
pD.lac <- filter(pD, Condition=="L")
cells <- as.character(pD.lac$barcode)
m.lac.all <- m[,cells]
keep <- rowMeans(m.lac.all) > 0.1
m.lac <- m.lac.all[keep,]
fD.lac <- fD[keep,]
#Normalization
clusters.lac <- quickCluster(m.lac)
pD.lac$sf <- computeSumFactors(m.lac,clusters=clusters.lac,positive=FALSE)
m.lac.norm <- t(t(m.lac)/pD.lac$sf)


fPCA.lac <- log2(t(m.lac.norm)+1)
fPCA.lac <- scale(fPCA.lac,scale=TRUE,center=TRUE)
set.seed(rnd_seed)
tsn.lac <- Rtsne(fPCA.lac,perplexity=25)

pD.lac$tSNE1 <- tsn.lac$Y[,1]
pD.lac$tSNE2 <- tsn.lac$Y[,2]

pD.lac <- arrange(pD.lac, PassAll) %>%
    rename(QCpass=PassAll) %>%
    mutate(QCpass=ifelse(QCpass,"Pass","Fail"))
tsnPlot <- ggplot(pD.lac, aes(x=tSNE1, y=tSNE2,color=QCpass)) +
    geom_point(alpha=0.8,size=3) +
    scale_color_manual(values=c("grey","black")) +
    theme(legend.title=element_blank())

cairo_pdf("../paper/figures/S2.pdf")
tsnPlot
dev.off()

# library(dynamicTreeCut)
# cls <- dynamicCluster(t(m.lac),ds=2,output="blub",minSize=10)
# pD.lac$cluster <- as.factor(cls$cluster)
# tsnPlot <- ggplot(pD.lac, aes(x=tSNE1, y=tSNE2,color=as.factor(cluster))) +
#     geom_point(alpha=0.8,size=3) +
    #     theme_bw() 
    # tsnPlot
    # 
    # 
    # out <- data.frame("C1"=numeric(nrow(m.lac.all)))
    # 
    # for (clust in levels(pD.lac$cluster)[-1]) {
    #     expr <- rowMeans(m.lac.all[,pD.lac$cluster==clust])
    #     colname <- paste0("C",clust)
    #     out[,colname] <- expr
    # }
    # remvd <- out
    # rownames(remvd) <- rownames(m.lac.all)
    # 
#     RELOAD DATA
# dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
# m <- dataList[[1]]
# pD <- dataList[[2]]
# fD <- dataList[[3]]
# relevel pD Condition to get it in logical sequence
# rnd_seed <- 300
# 
# pD$Condition <- mapvalues(pD$Condition, from=c("V","P","L","I"),
#                           to=c("NP", "G",
#                                "L", "PI"))
# pD$Condition <- factor(pD$Condition, levels=c("NP","G","L","PI"))
# pD$SampleID <- mapvalues(pD$SampleID, from=c("V1","V2","P1","P2","L1","L2","I1","I2"),
#                           to=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2"))
# pD$SampleID <- factor(pD$SampleID,levels=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2")) 
# 
# pD  <- filter(pD, PassAll)
# tach X
# 
# stopifnot(identical(rownames(remvd),rownames(m)))
# new.m <- as.matrix(cbind(m,remvd))
# add.pD <- pD[1:5,] %>% mutate(barcode=c("C1","C2","C3","C4","C5")) %>%
#     mutate(Condition="L",cluster=12)
# 
# new.pD <- rbind(pD,add.pD)
# 
# 
# A more or less fair representation of the lactation samples
# pD.lac <- filter(new.pD, Condition=="L")
# cells <- as.character(pD.lac$barcode)
# m.lac <- new.m[,cells]
# keep <- rowMeans(m.lac) > 0.1
# m.lac <- m.lac[keep,]
# fD.lac <- fD[keep,]
# Normalization
# clusters.lac <- quickCluster(m.lac)
# pD.lac$sf <- computeSumFactors(m.lac,clusters=clusters.lac,positive=FALSE)
# m.lac.norm <- t(t(m.lac)/pD.lac$sf)
# 
# 
# fPCA.lac <- log2(t(m.lac)+1)
# fPCA.lac <- scale(fPCA.lac,scale=TRUE,center=TRUE)
# set.seed(rnd_seed)
# tsn.lac <- Rtsne(fPCA.lac,perplexity=25)
# 
# pD.lac$tSNE1 <- tsn.lac$Y[,1]
# pD.lac$tSNE2 <- tsn.lac$Y[,2]
# 
# pD.lac <- arrange(pD.lac,cluster)
# tsnPlot <- ggplot(pD.lac, aes(x=tSNE1, y=tSNE2,color=cluster)) +
#     geom_point(alpha=0.8,size=3) +
# scale_color_manual(values=pal) +
    #     theme_bw() 
    # tsnPlot
    # 
    # library(dynamicTreeCut)
    # 
    # out <- data.frame("C1"=numeric(nrow(m.lac)))
    # 
    # for (clust in levels(pD.lac$cluster)[-1]) {
    #     expr <- rowMeans(m.lac[,pD.lac$cluster==clust])
    #     colname <- paste0("C",clust)
    #     out[,colname] <- expr
    # }
    # x <- out
    # rownames(x) <- rownames(m.lac)

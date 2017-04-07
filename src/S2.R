############################
#
# Figure S2
#
############################

library(scran)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Rtsne)
source("functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
rnd_seed <- 300

# Subset to lactation sample
pD.lac <- filter(pD, Condition=="L")
cells <- as.character(pD.lac$barcode)
m.lac.all <- m[,cells]
keep <- rowMeans(m.lac.all) > 0.1
m.lac <- m.lac.all[keep,]
fD.lac <- fD[keep,]

# Normalize
clusters.lac <- quickCluster(m.lac)
pD.lac$sf <- computeSumFactors(m.lac,clusters=clusters.lac,positive=FALSE)
m.lac.norm <- t(t(m.lac)/pD.lac$sf)

# Compute tSNE
fPCA.lac <- log2(t(m.lac.norm)+1)
fPCA.lac <- scale(fPCA.lac,scale=TRUE,center=TRUE)
set.seed(rnd_seed)
tsn.lac <- Rtsne(fPCA.lac,perplexity=25)

pD.lac$tSNE1 <- tsn.lac$Y[,1]
pD.lac$tSNE2 <- tsn.lac$Y[,2]

# Plot
pD.lac <- arrange(pD.lac, PassAll) %>%
    rename(QCpass=PassAll) %>%
    mutate(QCpass=ifelse(QCpass,"Pass","Fail"))
tsnPlot <- ggplot(pD.lac, aes(x=tSNE1, y=tSNE2,color=QCpass)) +
    geom_point(alpha=0.8,size=3) +
    scale_color_manual(values=c("grey","black")) +
    theme(legend.title=element_blank())

# cairo_pdf("../paper/figures/S2.pdf")
tsnPlot
# dev.off()


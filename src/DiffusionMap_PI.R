# Diffusion map predicted from PI on map from NP and G

library(destiny)
library(scran)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm_clustered_clean.rds")
dmList <- readRDS("../data/Robjects/DiffusionMap_Luminal.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
rm(dataList)

# Remove QC-fails,outlier and immune cells
excludeClustComb <- c("Bsl-G1","Bsl","Myo","Prc")
keep <- pD$keep & pD$Condition=="PI" & !(pD$SubCluster %in% excludeClustComb)
m.vp <- m[,keep]
pD.vp <- pD[keep,]

# Select genes that were used for diffusion map
keep <- dmList$genes
m.vp <- m.vp[keep,]

# Normalize 
m.vp <- t(t(m.vp)/pD.vp$sf)

# Prepare expression matrix
m.vp <- t(log(m.vp+1))

# Compute diffusion map
dm.pred <- dm_predict(dmList$dm,m.vp)
out <- data.frame(as.matrix(dm.pred[,1:4]),
		  barcode=pD.vp$barcode)
write.csv(out,"../data/Robjects/dm_luminal_PI.csv",row.names=FALSE)

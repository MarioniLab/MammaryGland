# Diffusion map computed on all cells from NP and G

library(destiny)
library(scran)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_QC_norm_clustered_clean.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
rm(dataList)

# Remove QC-fails,outlier and immune cells
m <- m[,pD$keep]
pD <- pD[pD$keep,]

# ---- BasalAndLuminalNPandG ----

condComb <- c("NP","G")
cells <- list("all","luminal")

for (cell in cells) {
    if(cell=="luminal") {
	excludeClustComb <- c("Bsl-G1","Bsl","Myo","Prc")
    } else{
	excludeClustComb <- NULL
    }

    keepCells <- pD$Condition %in% condComb & !(pD$SubCluster %in% excludeClustComb)
    m.vp <- m[,keepCells]
    pD.vp <- pD[keepCells,]
    pD.vp$SubCluster <- factor(pD.vp$SubCluster)

    # Remove genes below 0.1 mean
    keep <- rowMeans(m.vp)>0.01
    m.vp <- m.vp[keep,]

    # Normalize 
    m.vp <- t(t(m.vp)/pD.vp$sf)

    # Highly variable genes 
    var.des <- trendVar(log2(m.vp+1),trend="semiloess")
    var.out <- decomposeVar(log2(m.vp+1),var.des)
    o <- order(var.out$mean)
    hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >=0.5),]

    # Prepare expression matrix
    m.vp <- m.vp[rownames(hvg.out),]
    m.vp <- t(log(m.vp+1))

    # Compute diffusion map
    set.seed(rnd_seed)
    dm <- DiffusionMap(m.vp, n_eigs=20, rotate=TRUE)
    dms <- eigenvectors(dm)[,1:4]
    dms <- data.frame(dms,
		      barcode=pD.vp$barcode)

    #     library(dplyr)
    #     pD.vp <- left_join(pD.vp,dms,by="barcode")
    #     ggplot(pD.vp, aes(DC1,DC2,color=branch)) +
    #         geom_point() +
    #         theme_bw()
    
    # branching and dpt for luminal cells
    if(cell=="luminal") {

    # save for predicting PI later
    saveRDS(list(genes=colnames(m.vp),dm=dm),file="../data/Robjects/secondRun_2500/DiffusionMap_Luminal.rds")
    #Define tips as the cells at the corners of the triangluar shape
    t1 <- which.max(dms[,2])
    t2 <- which.min(dms[,1])
    t3 <- which.max(dms[,1])

    # Compute Pseudotime and branching
    set.seed(rnd_seed)
    dpt <- DPT(dm, branching=TRUE, tips=c(t1,t2,t3))
    root <- which(dpt@tips[,1])[2]
    rootdpt <- paste0("DPT",root)

    # Rename branches
    branch <- dpt@branch[,1]
    branch[is.na(branch)]  <- "Intermediate"
    branch[branch==1] <- "Root"
    branch[branch==2] <- "Secretory lineage"
    branch[branch==3] <- "Hormone-sensing lineage"

    #add branches and pseudotime to pD
    dms$branch <- factor(branch,levels=c("Root","Intermediate",
					"Secretory lineage",
					"Hormone-sensing lineage"))
    dms$dpt<- dpt[["dpt"]]
    }

    # save
    fileName <- sprintf("../data/Robjects/secondRun_2500/dm_%s.csv",cell)
    write.csv(dms, file=fileName, row.names=FALSE)
}


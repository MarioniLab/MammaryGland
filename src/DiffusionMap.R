# Diffusion map computed on all cells from NP and G

library(destiny)
source("functions.R")

# Load Data
rnd_seed <- 300
dataList <- readRDS("../data/Robjects/ExpressionList_Clustered.rds")
m <- dataList[[1]]
pD <- dataList[[2]]

# Remove QC-fails,outlier and immune cells
keepCells <- pD$PassAll & !pD$isImmuneCell & !pD$isOutlier 
m <- m[,keepCells]
pD <- pD[keepCells,]

# ---- BasalAndLuminalNPandG ----

condComb <- c("NP","G")
cells <- list("all","luminal")

for (cell in cells) {
    if(cell=="luminal") {
	excludeClustComb <- c(6,7,9)
    } else{
	excludeClustComb <- NULL
    }

    keepCells <- pD$Condition %in% condComb & !(pD$cluster %in% excludeClustComb)
    m.vp <- m[,keepCells]
    pD.vp <- pD[keepCells,]

    # Remove genes below 0.1 mean
    keep <- rowMeans(m.vp)>0.1
    m.vp <- m.vp[keep,]

    # Normalize 
    m.norm <- t(t(m.vp)/pD.vp$sf)

    # Compute HVGs
    brennecke <- BrenneckeHVG(m.norm,suppress.plot=TRUE)

    # Prepare expression matrix
    exps <- m.norm[brennecke,]
    exps <- t(log(exps+1))

    # Compute diffusion map
    set.seed(rnd_seed)
    dm <- DiffusionMap(exps, n_eigs=20, k=50, rotate=TRUE)
    dms <- eigenvectors(dm)[,1:3]
    dms <- data.frame(dms,
		      barcode=pD.vp$barcode)
    
    # branching and dpt for luminal cells
    if(cell=="luminal") {

    #Define tips as the cells at the corners of the triangluar shape
    t1 <- which.min(dms[,2])
    t2 <- which.min(dms[,1])
    t3 <- which.max(dms[,1])

    # Compute Pseudotime and branching
    set.seed(rnd_seed)
    dpt <- DPT(dm, branching=TRUE, tips=c(t1,t2,t3))
    root <- which(dpt@tips[,1])[3]
    rootdpt <- paste0("DPT",root)

    # Rename branches
    branch <- dpt@branch[,1]
    branch[is.na(branch)]  <- "Intermediate"
    branch[branch==1] <- "Root"
    branch[branch==2] <- "Hormone-sensing lineage"
    branch[branch==3] <- "Secretory lineage"

    #add branches and pseudotime to pD
    dms$branch <- factor(branch,levels=c("Root","Intermediate",
					"Secretory lineage",
					"Hormone-sensing lineage"))
    dms$dpt<- dpt[["dpt"]]
    }

    # save
    fileName <- sprintf("../data/Robjects/dm_%s.csv",cell)
    write.csv(dms, file=fileName, row.names=FALSE)
}


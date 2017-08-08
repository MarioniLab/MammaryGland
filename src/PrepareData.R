# Script to prepare cellranger data for downstream analysis

library(cellrangerRkit)
library(dplyr)

# ---- ReadData ----

# Read in output from cell ranger using cellrangerRkit
gene_bc_matrix <- load_cellranger_matrix("../data/CellRangerData/secondRun_2500/MammaryGland",
					 genome="mm10")
pDat <- data.frame(pData(gene_bc_matrix))
cDat <- as.matrix(exprs(gene_bc_matrix))

#reduce size of matrix
keep <- rowSums(cDat) > 1
cDat <- cDat[keep,]
fDat <- data.frame(fData(gene_bc_matrix))[rownames(cDat),]

# ---- Formatting ----

# Add more info to phenotype Data
pDat <- mutate(pDat, SeqID=substr(barcode,18,18)) %>%
	mutate(barcode=as.character(barcode)) %>%
        mutate(SeqID=mapvalues(SeqID,c("1","2","3","4","5","6","7","8"),
				 c("A7","B7","Myo","D7","E7","F7","G7","H7"))) %>%
        mutate(SampleID=mapvalues(SeqID,c("A7","B7","Myo","D7","E7","F7","G7","H7"),
				  c("PI1","PI2","L1","L2","G1","G2","NP1","NP2"))) %>%
	mutate(Replicate=as.factor(substr(SampleID,nchar(SampleID),nchar(SampleID))),
	       Condition=factor(substr(SampleID,1,nchar(SampleID)-1),levels=c("NP","G","L","PI"))) %>%
	mutate(SampleID=factor(SampleID,levels=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2"))) 

# Add more info to the feature Data
mitoGenes <- read.table("../data/miscData/MitoGenes.txt")
tfCheck <- read.table("../data/miscData/TFcheckpoint_WithENSID.tsv",
		header=TRUE, sep="\t")

fDat$Mitochondrial <- fDat$id %in% mitoGenes$V1
fDat$TranscriptionFactor <- fDat$id %in% tfCheck$ensembl_gene_id


# Save data
stopifnot(identical(rownames(fDat),rownames(cDat)) & identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../data/Robjects/secondRun_2500/ExpressionList.rds")

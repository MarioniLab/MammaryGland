#
# Script to prepare cellranger data for downstream analysis
#

library(cellrangerRkit)
library(dplyr)

# Read in output from cell ranger using cellrangerRkit
gene_bc_matrix <- load_cellranger_matrix("../data/CellRangerData/MammaryGland",
					 genome="mm10")
pDat <- data.frame(pData(gene_bc_matrix))
fDat <- data.frame(fData(gene_bc_matrix))
cDat <- as.matrix(exprs(gene_bc_matrix))

# Add more info to phenotype Data
pDat <- mutate(pDat, SeqID=substr(barcode,18,18)) %>%
	       mutate(SeqID=mapvalues(SeqID,c("1","2","3","4","5","6","7","8"),
					 c("A1","B1","C1","D1","E1","F1","G1","H1"))) %>%
        mutate(SampleID=mapvalues(SeqID,c("A1","B1","C1","D1","E1","F1","G1","H1"),
				  c("NP1","G2","PI2","G1","L2","L1","NP2","PI1"))) %>%
	mutate(Replicate=as.factor(substr(SampleID,nchar(SampleID),nchar(SampleID))),
	       Condition=factor(substr(SampleID,1,nchar(SampleID)-1),levels=c("NP","G","L","PI"))) %>%
	mutate(SampleID=factor(SampleID,levels=c("NP1","NP2","G1","G2","L1","L2","PI1","PI2"))) 

# Add more info to the feature Data
mitoGenes <- read.table("../data/miscData/MitoGenes.txt")
tfCheck <- read.table("../data/miscData/TFcheckpoint_WithENSID.tsv",
		header=TRUE, sep="\t")

fDat <- mutate(fDat,
	       Mitochondrial=id %in% mitoGenes$x,
	       TranscriptionFactor=id %in% tfCheck$ensembl_gene_id)


# Save data
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../data/Robjects/ExpressionList.rds")

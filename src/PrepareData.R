library(cellrangerRkit)
gene_bc_matrix <- load_cellranger_matrix("../data/CellRangerData/MammaryGland/", genome="mm10")
pDat <- data.frame(pData(gene_bc_matrix))
fDat <- data.frame(fData(gene_bc_matrix))
cDat <- as.matrix(exprs(gene_bc_matrix))

#Add more info to phenotype Data
library(dplyr)
pDat <- mutate(pDat, SeqID=substr(barcode,18,18)) %>%
	       mutate(SeqID=mapvalues(SeqID,c("1","2","3","4","5","6","7","8"),
					 c("A1","B1","C1","D1","E1","F1","G1","H1"))) %>%
        mutate(SampleID=mapvalues(SeqID,c("A1","B1","C1","D1","E1","F1","G1","H1"),
				  c("V1","P2","I2","P1","L2","L1","V2","I1"))) %>%
	mutate(Replicate=as.factor(substr(SampleID,2,2)),
	       Condition=as.factor(substr(SampleID,1,1)))

#Add more info to the feature Data
mitoGenes <- read.table("../data/miscData/MitoGenes.txt")
fDat <- mutate(fDat,
	       Mitochondrial=id %in% mitoGenes$x)


DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../data/Robjects/ExpressionList.rds")

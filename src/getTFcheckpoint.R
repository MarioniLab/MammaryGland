# Scritp to download and format TFcheckpoint database

library(dplyr)
library(biomaRt)

# ---- Download ----
download.file("http://www.tfcheckpoint.org/data/TFCheckpoint_download_180515.txt", destfile="../data/miscData/TFcheckpoint.txt")

tfcheckpoint <- read.csv("../data/miscData/TFcheckpoint.txt", sep="\t", row.names=NULL, stringsAsFactors=FALSE,
			   header=TRUE)

# Keep only tfs with human entrezID & literature support
tfcheckpoint <- tfcheckpoint[tfcheckpoint$entrez_mouse!=0,]
tfcheckpoint <- tfcheckpoint[tfcheckpoint$DbTF!=0,]

# ---- formating ----
# Add column with ENS ID
ensembl <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")

filt <- tfcheckpoint$entrez_mouse
res <- getBM(mart=ensembl,
	     filters="entrezgene",
	     values=filt,
	     attributes=c("entrezgene","ensembl_gene_id"))


colnames(res) <- c("entrez_mouse","ensembl_gene_id")
tfcheckpoint.withENSG <- full_join(res,tfcheckpoint,by="entrez_mouse")

# Save
write.table(tfcheckpoint.withENSG,"../data/miscData/TFcheckpoint_WithENSID.tsv",row.names=FALSE,quote=FALSE,
	  sep="\t")


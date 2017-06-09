#!/bin/bash

# download gtf file from ensembl
wget -O ../data/miscData/mm10.gtf.gz "ftp://ftp.ensembl.org/pub/release-89/gtf/mus_musculus/Mus_musculus.GRCm38.89.chr.gtf.gz" 
gunzip ../data/miscData/mm10.gtf.gz

# extract mitochondrial genes
grep ^MT ../data/miscData/mm10.gtf | grep -o 'ENSMUSG[0-9]*' | uniq > ../data/miscData/MitoGenes.txt

# delete gtf
rm ../data/miscData/mm10.gtf 

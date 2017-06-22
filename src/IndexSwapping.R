library(ggplot2)
library(cowplot)

# Computation of p-values based on hypergeometric test
compare <- function(barcodes, samples) {
    out <- NULL
    ids <- levels(samples)
    combs <- combn(ids,m=2, simplify=FALSE)
    for (i in seq_along(combs)) {
	comb <- combs[[i]]
	s1 <- comb[1]
	s2 <- comb[2]
	bc1 <- as.character(barcodes[samples==s1])
	bc2 <- as.character(barcodes[samples==s2])
	x <- length(intersect(bc1,bc2))
	m <- length(bc1)
	n <- 750000
	k <- length(bc2)
	p.val <- phyper(q=x-1,m=m,n=n,k=k,lower.tail=FALSE)
	tmp <- data.frame(s1=s1,
			  s2=s2,
			  n1=m,
			  n2=k,
			  shared=x,
			  p.val=p.val)
	out <- rbind(out,tmp)
    }
    return(out)
}

# Load Data
dList <- readRDS("../data/Robjects/secondRun_2500/ExpressionList_Clustered.rds")
pD <- dList[[2]]

# Extract barcodes by removing sampleID from bc
pD <- pD[pD$PassAll & !pD$isImmuneCell & !pD$isOutlier,]
barcodes <- as.character(pD$barcode)
spt <- strsplit(barcodes, split = "-", fixed = T)
pD$sample <- sapply(spt, function(x) x[2])
pD$bcs <- sapply(spt, function(x) x[1])
pD.add <- data.frame(bcs=names(table(pD$bcs)),
		     bc.obs=(table(pD$bcs)))
pD.add <- pD.add[,-2]
pD <- dplyr::left_join(pD,pD.add)

# ------- More Barcodes are shared than expected -------------

# P1
p1 <- ggplot(pD, aes(x=SampleID, fill=as.factor(bc.obs.Freq))) +
    geom_bar() +
    #     ggtitle("Samples contain many shared barcodes") +
    scale_fill_discrete(name="Times barcode observed") +
    theme(legend.position="bottom",
	  legend.direction="horizontal")

ggsave("../submissions/resubmission/img/p1.pdf",p1)

# P2
compDf <- compare(pD$bcs, pD$SampleID)
p2 <- ggplot(compDf, aes(x=shared, y=-log10(p.val))) +
    geom_point() +
    xlab("# shared Barcodes") +
    ylab("-log10(P)") +
    geom_hline(yintercept=2,lty="dashed",color="red") 
    #     ggtitle("Samples share more barcodes than expected by chance") 

ggsave("../submissions/resubmission/img/p2.pdf",p2)

bcTab <- as.matrix(table(pD$bcs,pD$cluster))
dupBcs <- pD[pD$bc.obs.Freq>1,"bcs"]
bcTab <- bcTab[pD[pD$bcs %in% dupBcs,"bcs"],]
freq <-  apply(bcTab,1,function(x) ifelse(max(x)>1,max(x),0))/apply(bcTab,1,sum) * 100
freq <- freq[unique(names(freq))]
pD.add <- data.frame("bcs"=names(freq),
		   "freqSameCluster"=freq)
pD <- dplyr::left_join(pD,pD.add,by="bcs")

p3 <- ggplot(pD, aes(x=freqSameCluster)) +
    geom_bar(color="black",fill="grey") +
    xlab("Fraction of cells per shared barcode in the same cluster") +
    ylab("Count") 
    #     ggtitle("Cells that share barcodes cluster together")

ggsave("../submissions/resubmission/img/p3.pdf",p3)

library(dplyr)

ratioDf <- group_by(pD,bcs) %>%
    mutate(largestCell=UmiSums==max(UmiSums)) %>%
    mutate(largestSample=SampleID[largestCell]) %>%
    mutate(ratio=max(UmiSums)/((sum(UmiSums)-max(UmiSums))/(n()-1))) %>%
    ungroup() %>%
    filter(largestCell)

p4 <- ggplot(ratioDf, aes(x=ratio, fill=largestSample)) +
    geom_density() +
    xlab("Library size ratio of largest cell versus mean of remaining cells") 
    #     ggtitle("Barcode groups contain one large cell and many small ones")

ggsave("../submissions/resubmission/img/p4.pdf",p4)

# ---------- Lactation sample did not contain real cells --------------

p5 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=Condition)) +
    theme_void(base_size=12) +
    geom_point(size=1.5) +
    geom_line(aes(group=bcs, color=NULL), lty="dashed",size=0.3, alpha=0.7) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="bottom",legend.direction="horizontal",
	  legend.title=element_blank()) 

p5    
ggsave("../submissions/resubmission/img/p5.pdf",p5)

fPlot <- filter(pD, bc.obs.Freq > 1) %>%
    group_by(bcs) %>%
    summarise(potOrigin=Condition[which.max(UmiSums)],
	      targets=paste(sort(SampleID[-which.max(UmiSums)]),collapse="-")) %>% 
    mutate(pattern=paste(potOrigin,targets,"_"))


p6 <- ggplot(fPlot, aes(x=targets, fill=potOrigin)) +
    geom_bar() +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    xlab("Sample ID of smaller cells") +
    scale_fill_discrete(name="Largest cell from") 

ggsave("../submissions/resubmission/img/p6.pdf",p6)

test <- filter(pD, bc.obs.Freq > 1) %>%
    group_by(bcs) %>%
    summarise(conts=any(Condition=="L"),
	      LacSmall=Condition[which.max(UmiSums)]!="L")

trueL <- group_by(pD, bcs) %>%
    summarise(conts=Condition[which.max(UmiSums)]) %>%
    filter(conts=="L") %>%
    .$bcs

sbst <- filter(pD, bcs %in% trueL) %>%
    mutate(cluster=paste0("C",cluster))

p7 <- ggplot(sbst, aes(x=bcs,y=UmiSums /1000,fill=SampleID)) +
    geom_bar(stat="identity",position="dodge") +
    facet_grid(~cluster, scales="free") +
    scale_color_brewer("Paired") +
    ylab("Library Size [x1000]") +
    xlab("Barcode") +
    coord_flip()

ggsave("../submissions/resubmission/img/p7.pdf",p7)



# Load Data
sbst <- pD[pD$Condition !="L" & !(pD$bcs %in% trueL), !(colnames(pD) %in% c("sample","bcs","bc.obs.Freq"))]
barcodes <- as.character(sbst$barcode)
spt <- strsplit(barcodes, split = "-", fixed = T)
sbst$sample <- sapply(spt, function(x) x[2])
sbst$bcs <- sapply(spt, function(x) x[1])
sbst.add <- data.frame(bcs=names(table(sbst$bcs)),
		     bc.obs=(table(sbst$bcs)))
sbst.add <- sbst.add[,-2]
sbst <- dplyr::left_join(sbst,sbst.add)


p8 <- ggplot(sbst, aes(x=SampleID, fill=as.factor(bc.obs.Freq))) +
    geom_bar() +
    #     ggtitle("Samples contain many shared barcodes") +
    scale_fill_discrete(name="Times barcode observed") +
    theme(legend.position="bottom",
	  legend.direction="horizontal")
compDf <- compare(sbst$bcs, sbst$SampleID)
p9 <- ggplot(compDf, aes(x=shared, y=-log10(p.val))) +
    geom_point() +
    xlab("# shared Barcodes") +
    ylab("-log10(P)") +
    geom_hline(yintercept=2,lty="dashed",color="red") 
    #     ggtitle("Samples share more barcodes than expected by chance")

p10 <- plot_grid(p8,p9,nrow=2)

ggsave("../submissions/resubmission/img/p10.pdf",p10)

# end



























## Amount of reads that bled into the lactation sample

So far it is seems fairly convincing that the lactation samples almost inclusively possess bleed-throughs from other samples. In theory we expect a large cell with many reads artificially causing one or two smaller cells with the same barcode in both lactation samples (and maybe other samples). We can now get an overview of the extent of index swapping in terms of the fraction of total reads of the biggest cell versus second biggest cell. In the following plot only barcodes that come up in at least one lactation sample are selected. The plot shows the distributions of the aforementioned fractions for these barcodes. Note that this overestimates the fraction of reads that bled into the lactation sample as some of the reads are likely attributable to a background of mRNA floating around in the sample. A preciser (and more time consuming) way of estimating the amount of index swapping would be to look at the UMIs. However, the rough order of magnitutex appears to be similar to what has been reported.

```{r echo=FALSE, message=FALSE}
frcts <- sapply(bleedBcs, function(x) {
		   sbst <- pD[pD$bcs==x,]
		   frct <- sort(sbst$UmiSums, decreasing=TRUE)[2]/max(sbst$UmiSums)
	  })

conds <- sapply(bleedBcs, function(x) {
		   sbst <- pD[pD$bcs==x,]
		   sbst.max <- sbst[sbst$UmiSums==max(sbst$UmiSums),"Condition"]
	  })
dat <- data.frame("Fraction"=frcts,
		  "Condition"=conds)
dat <- dat[dat$Condition!="L",]

ggplot(dat, aes(x=Fraction)) +
       geom_histogram() +
       theme_bw() +
       ggtitle("Library size fraction of biggest cell compared to second biggest cell") +
       xlim(c(0,1))
```

# What's the impact on our biological conclusions

## Results form F1

The overall clustering will be affected to some extent, as a substantial amount of cells cannot be trusted and would need to be removed.
If we tabulate cluster assignment (rows) versus the frequency a barcode appears in the dataset (cols) it appears that most clusters contain a majority of cells with unique barcodes which is reassuring. 

```{r echo=FALSE, message=FALSE}
pD$cluster <- factor(pD$cluster)
table(pD$cluster,pD$bc.obs.Freq)
```

__However, cluster 7 and cluster 8 only contain barcodes that appear in multiple samples.__
My feeling is that these clusters stem from a few cells that were actually present in the lactation sample and whose cDNA was then highly present in the sequencing (due to equal mixing of all cDNA libraries) and therefore was more prone to index swapping.
In fact, for all of these groups of barcodes the cell with the largest library size came from the Lactation samples, suggesting index swapping from lactation to the other samples.
```{r echo=FALSE, message=FALSE}
library(RColorBrewer)
bleedBcs <- unique(pD[pD$cluster %in% c(7,8),"bcs"])
sbst <- pD[pD$bcs %in% bleedBcs,]
ggplot(sbst, aes(x=bcs,y=UmiSums,fill=SampleID)) +
    geom_bar(stat="identity",position="dodge") +
    facet_grid(~cluster, scales="free") +
    scale_color_brewer("Paired") +
    ylab("Library Size") +
    xlab("Barcode") +
    coord_flip()
```
Conclusively, if these clusters represent real cell types they only appear in the lactation sample. I think that cluster 7 might genuinely be a myoepithelial cell, the expression profile is quite convincing. Yet, obviously these clusters now represent only few 'true' cells and we should be very skeptical about them.

## Results from F2 and F3

I don't think that our main biological conclusion will change much, apart from the fact that I am not sure any more whether cluster 8 is actually real. In anyway if it was it stems from the lactation sample and would thus not be included in the trajectory. Otherwise the conclusion that the luminal compartment has one common progenitor is not going to change. Also differential expression along the trajectories will be more or less fine, apart from effects that stem from C8.

## Results from F4

Here again the main conclusion remains the same, that is C4 is a post-parity specific progenitor cell. The only conclusion that would change is F4d, which is that the expression of milk genes is equal between C4-L and C4-PI, as C4-L is not real.

# Conclusions

The impact of the index swapping issue of the ExAmp chemistry from the HiSeq4000 appears to be very severe for our 10x dataset. Having said that, the severity appears to stem from the two lactation samples.

From my current understanding both lactation samples were essentially devoid of cells (or cells didn't lyse).
The only reads that truly come from this sample are coming from 'background' mRNA that was floating in the medium from dead/bursted cells.
This leads to assignment of reads equally to all barcodes.
On top of that other samples bled through by index swapping.
This leads to cellular barcodes in the lactation sample containing reads from some other cell and background mRNA.
For the CellRanger pipeline it is then essentially impossible to distinguish these barcodes from genuine cells, especially as this sample did (almost) not contain any. 
The other samples instead contained enough genuine cells to allow the distinction from barcodes that only contained background + index swapping reads.
Hence, if the lactation samples are removed the number of non-unique barcodes per sample drops almost to zero.
Nonetheless, some cells are most likely still an artifact, especially cells in C7 and C8.

```{r echo=FALSE, message=FALSE, warning=FALSE}
pD <- pD[pD$PassAll & !pD$isImmuneCell & !pD$isOutlier & pD$Condition!="L",-c(25,26)]
pD.add <- data.frame(bcs=names(table(pD$bcs)),
		     bc.obs=(table(pD$bcs)))
pD.add <- pD.add[,-2]
pD <- left_join(pD,pD.add)

ggplot(pD, aes(x=SampleID, fill=as.factor(bc.obs.Freq))) +
    geom_bar() +
    theme_bw() +
    ggtitle("Number of unique barcodes omitting the Lactation sample")+
    scale_fill_discrete(name="# Barcode observed")
```


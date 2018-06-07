load("data/data.RData")
# dataList <- readRDS("data.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]


# Branch Data
out <- readRDS("data/BranchDEList.rds")

# both matrices
m.hrm <- out[[1]][["mSmooth"]]
m.alv <- out[[2]][["mSmooth"]]

# Results
res.hrm <- out[[1]][["Results"]]
res.alv <- out[[2]][["Results"]]

top.same <- as.character(out[["genes.sameGrad"]])[1:200]
top.diff <- as.character(out[["genes.diffGrad"]])[1:200]

plotGeneDist <- function(m, pD, fD, genes, colorBy="Condition") {
    # function to plot gene expression as boxplots
    require(reshape)
    feat <- fD[fD$symbol %in% genes,"id"]
    exps <- data.frame(as.vector(t(m[feat,])))
    colnames(exps) <- feat
    exps$barcode <- colnames(m)
    pltDat <- left_join(select_(pD, "barcode", colorBy),exps) %>% melt(.,id=c("barcode",colorBy))
    pltDat[,colorBy] <- as.factor(pltDat[,colorBy])
    plt <- ggplot(pltDat, aes_string(x=colorBy,y="value",fill=colorBy)) +
	geom_boxplot() +
	#         geom_point(position="jitter",alpha=0.2,shape=19) +
	ggtitle(genes) +
	scale_fill_manual(values=levels(pD$Colors)) +
	guides(fill=FALSE) 
    return(plt)
}

ggname <- function(x) {
    if (class(x) != "character") {
        return(x)
    }
    y <- sapply(x, function(s) {
        if (!grepl("^`", s)) {
            s <- paste("`", s, sep="", collapse="")
        }
        if (!grepl("`$", s)) {
            s <- paste(s, "`", sep="", collapse="")
        }
    }
    )
    y 
}

plotPseudotimeTrend <- function(input) {
    if (input$list=="All") {
	feature <- input$Gene2
    } else
	if (input$list=="Top genes with similar gradient") {
	feature <- input$Gene3
    } else {
	feature <- input$Gene4
    }
    #extract expression for first branch
    p1 <- out[[1]][["pD"]] %>% arrange(DPTRank)
    fplot1 <- data.frame(t(m.hrm)[,feature],
			check.names=FALSE)
    colnames(fplot1) <- feature
    fplot1$barcode <- as.character(p1$barcode)
    fplot1$raw <- as.vector((out[[1]][["m"]][feature,fplot1$barcode]))
    colnames(fplot1)[ncol(fplot1)] <- paste0("raw",feature)
    fplot1 <- join(fplot1,p1,by="barcode") %>%
	mutate(dptNorm=dpt/max(dpt))

    #extract expression for second branch
    p2 <- out[[2]][["pD"]] %>% arrange(DPTRank)
    fplot2 <- data.frame(t(m.alv)[,feature],
			check.names=FALSE)
    colnames(fplot2) <- feature
    fplot2$barcode <- as.character(p2$barcode)
    fplot2$raw <- as.vector((out[[2]][["m"]][feature,fplot2$barcode]))
    colnames(fplot2)[ncol(fplot2)] <- paste0("raw",feature)
    fplot2 <- join(fplot2,p2,by="barcode") %>%
	mutate(dptNorm=dpt/max(dpt))

    #set colorscale
    clustCol <- levels(fplot2$Colors)[levels(fplot2$SubCluster)[c(-12,-17,-18,-19,-20)] %in% unique(c(as.character(fplot2$SubCluster),as.character(fplot1$SubCluster)))]

    fplot3 <- rbind(fplot1,fplot2)
    #Plot dpt-dependent expression for feature
    pnts <- ggname(paste0("raw",feature))
    lns <- ggname(feature)
    p <- ggplot() +
	geom_point(size=0.8,data=fplot3,aes_string(x="dptNorm",y=pnts, color="SubCluster")) +
	geom_line(data=fplot1,aes_string(x="dptNorm",y=lns),lty="dashed") +
	geom_line(data=fplot2,aes_string(x="dptNorm",y=lns)) +
	ggtitle(feature) +
	ylab("Log-Expression") +
	scale_color_manual(values=clustCol) +
	scale_x_continuous(breaks=c(0,1)) +
	theme(legend.position="bottom",
	      legend.direction="horizontal",
	      legend.title=element_blank(),
	      plot.title=element_text(face="italic")) +
	guides(colour = FALSE) +
	xlab("Pseudotime")

    # create a dummy legend for dasehd/solid line
    dummyD <- data.frame(x=rnorm(10),y=rnorm(10),
			 branch=c(rep("Secretory lineage",5),rep("Hormone-sensing lineage",5)))
    dummyP2 <- ggplot(dummyD,aes(x,y,lty=branch)) +
	geom_line() +
	theme(legend.position="bottom",
	      #           legend.direction="horizontal",
	      legend.title=element_blank()) +
	scale_linetype_manual(values=c("dashed","solid"))
    legs1 <- get_legend(dummyP2)

    plot_grid(p,legs1,rel_heights=c(1,0.1),nrow=2)
}

plotTriangleExpression <- function(input) {
    if (input$list=="All") {
	feature <- input$Gene2
    } else
	if (input$list=="Top genes with similar gradient") {
	feature <- input$Gene3
    } else {
	feature <- input$Gene4
    }

    # Genes to plot for trends
    exps <- m[fD$symbol==feature,]


    # setup df
    fPlot <- data.frame(exps,
			barcode=colnames(m))
    pD.sub <- pD[!is.na(pD$DC1),]
    fPlot <- right_join(fPlot, pD.sub[,c("barcode","DC1","DC2")], by="barcode")

    #     browser()
    fPl <- arrange(fPlot, exps)
    p <- ggplot(fPl, aes(x=DC1,y=DC2, color=exps)) +
	geom_point(size=4, pch=20) +
	scale_color_viridis() +
	ggtitle(feature) +
	theme_void() +
	theme(legend.position="bottom",
	      legend.direction="horizontal",
	      legend.title=element_blank(),
	      plot.title=element_text(face="italic"))

     return(p)
}

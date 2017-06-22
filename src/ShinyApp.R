library(topGO)
library(org.Mm.eg.db)
library(shiny)
library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
library(DT)
source(file.path("functions.R"))

dataList <- readRDS(file.path("../data/Robjects/secondRun_2500/ExpressionList_Clustered.rds"))
deList <- readRDS(file.path("../data/Robjects/secondRun_2500/DEList.rds"))
#Commence
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- pD$PassAll & ! (pD$isImmuneCell | pD$isOutlier)
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep <- rowSums(m)>10
m <- m[keep,]
fD <- fD[keep,]

#Normalization
m <- t(t(m)/pD$sf)

#Add cell surface annotation
library(scLVM)
surface <- getEnsembl("GO:0009986")
cellCycle <- getEnsembl("GO:0007049")
tfs <- filter(fD, TranscriptionFactor) %>% .$id

deAdd <- select(fD, id,symbol) 
colnames(deAdd) <- c("Gene","GeneSymbol")
deList <- lapply(deList, function(x) {
		      nw <- full_join(x,deAdd)
		      nw$Surface <- nw$Gene %in% surface
		      nw$CC <- nw$Gene %in% cellCycle
		      nw$TF <- nw$Gene %in% tfs
		      nw <- nw[,-2]
		      cols <- c("GeneSymbol",grep("top",colnames(nw),value=TRUE, #reorder cols
						  ignore.case=TRUE))
		      nw <- nw[,c(cols,setdiff(colnames(nw), cols))]})
names(deList) <- c(1,2,3,4,5,6,7,9,10)

ui <- shinyUI(fluidPage(

  titlePanel("Tabsets"),

  sidebarLayout(
    
    sidebarPanel(
    inputPanel(
	       selectInput("colorBy", label = "Color By",
			   choices= c("cluster","Condition"), selected="cluster"),
	       selectInput("cNumber", label = "Cluster",
			   choices= as.character(c(1,2,3,4,5,6,7,9,10)), selected="1"),
	       selectInput("ontology", label = "Ontology",
			   choices= c("BP","CC","MF"), selected="MF"),
	       radioButtons("order", label = "Order DE-Table",
			   choices= c("Up", "Down", "none"), selected="none"),
	       selectInput("against", label = "Compare Against",
			   choices= as.character(c(1,2,3,4,5,6,7,9,10)),
			   multiple=TRUE)
			   
	  ),plotOutput("distPlot1"),
	    plotOutput("distPlot2")
	

    ),
  
    mainPanel(
        tabsetPanel(
		    tabPanel("Plots",splitLayout(plotOutput("tSNE"), plotOutput("tSNE2")),
		dataTableOutput("table")),
		    tabPanel("Table",verbatimTextOutput("table2")),
		    tabPanel("GO-Analysis",tableOutput("goTable")))
      )
    )
  )
)

server <- function(input, output) {

dataSet <- reactive({
        out <- deList[[input$cNumber]]
	cols <- grepl(".vs.",colnames(out))
	out[,cols] <- round(out[,cols],digits=1)
    if(input$order=="none") {
	if(length(input$against)!=0) {
	cols <- paste0("logFC.vs.",input$against)
	colsP <- paste0("top.vs.",input$against)
	sbst <- out[,cols]
	ordr <- order(rowMax(as.matrix(cbind(out[,colsP],out[,colsP]))))
	out <- out[ordr,]
	}
	return(out)
    } else {
	cols <- paste0("logFC.vs.",input$against)
	colsP <- paste0("top.vs.",input$against)
	sbst <- out[,cols]
	if(input$order=="Up") {
	    sbst <- cbind(sbst,sbst) # prevent apply to crash for single cluster comparison
	    allUp <- apply(sbst, 1, function(x) sum(x>1)==length(x))
	    marker.up <- out[allUp,]
	    ordr <- order(rowMax(as.matrix(cbind(marker.up[,colsP],marker.up[,colsP]))))
	    marker.up <- marker.up[ordr,]
	    return(marker.up) }
	else {
	    sbst <- cbind(sbst,sbst) # prevent apply to crash for single cluster comparison
	    allDown <- apply(sbst, 1, function(x) sum(x<1)==length(x))
	    marker.down <- out[allDown,]
	    ordr <- order(rowMax(as.matrix(cbind(marker.down[,colsP],marker.down[,colsP]))))
	    marker.down <- marker.down[ordr,]
	    return(marker.down) }
    }})

output$tSNE <- renderPlot({
	tsnPlot <- ggplot(pD, aes_string(x="tSNE1", y="tSNE2", color=input$colorBy, shape="Replicate")) +
	    geom_point(size=1.5) +
	    scale_color_brewer(palette="Paired") +
	    theme_void()
	tsnPlot
    })

output$tSNE2 <- renderPlot({
	library(viridis)
	s <- input$table_rows_selected
	dat <- dataSet()
	gene <- dat[s,"GeneSymbol"]
	gn <- filter(fD, symbol %in% gene) %>% .$id
	expr <- log2(t(m)[,gn]+1)
	pD[,gn] <- expr
	pD <- arrange(pD, (pD[,gn]))
	tsnPlot <- ggplot(pD, aes_string(x="tSNE1", y="tSNE2", color=gn, shape="Replicate")) +
	    geom_point(size=1.5) +
	    scale_color_viridis() +
	    theme_void()
	tsnPlot
    })

output$distPlot1 <- renderPlot({
	s <- input$table_rows_selected
	dat <- dataSet()
	gene <- dat[s,"GeneSymbol"]
	     plotGeneDist(m, pD, fD, gene, colorBy="cluster")
       })

output$distPlot2 <- renderPlot({
	s <- input$table_rows_selected
	dat <- dataSet()
	pD <- filter(pD, cluster==input$cNumber)
	m <- m[,as.character(pD$barcode)]
	gene <- dat[s,"GeneSymbol"]
	     plotGeneDist(m, pD, fD, gene, colorBy="Condition")
       })

output$table <- renderDataTable({
        out <- dataSet()
	datatable(out, filter="top", selection=list(mode="single",
						    selected=1))
       })

output$table2 <- renderPrint({
     table(pD$cluster,pD$Condition)
       })

output$goTable <- renderTable({
    univrs <- fD$symbol
    out <- dataSet()
    top100 <- out$Gene[1:100]
    alG <- factor(as.numeric(univrs %in% top100))
    names(alG) <- univrs

    # ---- GOanalysis ----

    ## prepare Data for topGO
    ont <- input$ontology
    GO.data <- new("topGOdata", description="Lib GO",ontology=ont, allGenes=alG, 
		   annot=annFUN.org, mapping="org.Mm.eg.db",
		   nodeSize=10, ID="symbol")
    result.classic <- runTest(GO.data, algorithm="classic", statistic="Fisher")
    output <- GenTable(GO.data, Fisher.classic=result.classic, orderBy="Fisher.classic", topNodes=50, numChar=10000)
    return(output)
       })
}

shinyApp(ui = ui, server = server)

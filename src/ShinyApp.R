library(shiny)
library(dplyr)
library(gtools)
library(knitr)
library(ggplot2)
library(DT)
source(file.path("functions.R"))

rnd_seed <- 300
dataList <- readRDS(file.path("../data/Robjects/ExpressionList_Clustered.rds"))
deList <- readRDS(file.path("../data/Robjects/DEList.rds"))
#Commence
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]

#Pre-Filtering before DE-Analysis
# Cells
keepCells <- pD$PassAll & pD$cluster !=0
m <- m[,keepCells]
pD <- pD[keepCells,]

# Genes
keep1 <- rowSums(m!=0)>10
keep2 <- rowSums(m)>50
keep <- keep1 & keep2
m <- m[keep,]
fD <- fD[keep,]

#Normalization
m <- t(t(m)/pD$sf)

deAdd <- select(fD, id,symbol) 
colnames(deAdd) <- c("Gene","GeneSymbol")
deList <- lapply(deList, function(x) {
		      nw <- left_join(x,deAdd)
		      nw <- nw[,-2]})

ui <- shinyUI(fluidPage(

  titlePanel("Tabsets"),

  sidebarLayout(
    
    sidebarPanel(
    inputPanel(
	       selectInput("colorBy", label = "Color By",
			   choices= c("cluster","Condition"), selected="cluster"),
	       selectInput("cNumber", label = "Cluster",
			   choices= as.character(c(1,2,3,4,5,6,7,8,9,10)), selected="1"),
	       radioButtons("order", label = "Order DE-Table",
			   choices= c("Up against all", "Down against all", "none"), selected="none")
			   
	  ),plotOutput("distPlot1"),
	    plotOutput("distPlot2")
	

    ),
  
    mainPanel(
        tabsetPanel(
		    tabPanel("Plots",splitLayout(plotOutput("tSNE"), plotOutput("tSNE2")),
		dataTableOutput("table")),
		    tabPanel("Table",verbatimTextOutput("table2")))
      )
    )
  )
)

server <- function(input, output) {

dataSet <- reactive({
        out <- deList[[input$cNumber]]
	out[,2:10] <- round(out[,2:10],digits=2)
    if(input$order=="none") {
	return(out)
    } else {
	if(input$order=="Up against all") {
	    sbst <- out[,2:10] 
	    allUp <- apply(sbst, 1, function(x) sum(x>0)==length(x))
	    marker.up <- out[allUp,]
	    marker.up <- marker.up[order(apply(marker.up[,2:10],1,min),decreasing=TRUE),]
	    return(marker.up) }
	else {
	    sbst <- out[,2:10] 
	    allDown <- apply(sbst, 1, function(x) sum(x<0)==length(x))
	    marker.down <- out[allDown,]
	    marker.down <- marker.down[order(apply(marker.down[,2:10],1,max),decreasing=FALSE),]
	    return(marker.down) }
    }})

output$tSNE <- renderPlot({
	tsnPlot <- ggplot(pD, aes_string(x="tSNE1", y="tSNE2", color=input$colorBy, shape="Replicate")) +
	    geom_point(size=1.5) +
	    scale_color_brewer(palette="Paired") +
	    theme_bw()
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
	    theme_bw()
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
}

shinyApp(ui = ui, server = server)

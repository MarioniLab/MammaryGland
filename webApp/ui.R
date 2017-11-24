library(shiny)
library(plyr)
library(dplyr)
library(ggplot2)
library(gtools)
library(cowplot)
library(reshape)
source("helper.R")

shinyUI(fluidPage(

  navbarPage("Mammary Gland Development",
  tabPanel("About",
	   includeMarkdown("README.md")
	   ),
  tabPanel("Clusters",
  sidebarLayout(

    sidebarPanel(
	       selectInput("colorBy", label = "Color By",
			   choices= c("Cluster","InferredCellType","Condition"), selected="SuperCluster"),
	       selectInput("Gene", label = "Gene",
			   choices=mixedsort(fD$symbol), selected ="Csn2",
			   selectize=FALSE),
	      img(src="dendro.png",width="484px",height="363px")
    ),
  
    mainPanel(h4("Cluster"),
	     plotOutput("tSNE",width="25%",height="200px"),
	     h4("Gene Expression"),
	     splitLayout(plotOutput("tSNE2"),plotOutput("distPlot1")),
	     fluidRow(column(3,downloadButton('downloadPlot1')),
		      column(3, offset=3, downloadButton('downloadPlot2')))
      )
    )
  ),
 tabPanel("Differentiation trajectory",
	     sidebarLayout(
	       sidebarPanel(
			    selectInput("list", label= "Gene List", choices= c("All","Top genes with similar gradient","Top branch-specific")),
			    conditionalPanel("input.list=='All'",
			    selectInput("Gene2", label = "Gene",
					choices=mixedsort(rownames(m.hrm)), selected ="Csn2",
					selectize=FALSE)),
			    conditionalPanel("input.list=='Top genes with similar gradient'",
			    selectInput("Gene3", label = "Gene",
					choices=top.same,# selected ="Csn2",
					selectize=FALSE)),
			    conditionalPanel("input.list=='Top branch-specific'",
			    selectInput("Gene4", label = "Gene",
					choices=top.diff,# selected ="Csn2",
					selectize=FALSE)),
	      img(src="dendro.png",width="484px",height="363px")
	       ),

	       mainPanel(plotOutput("DiffMap",width="25%",height="200px"),
			 splitLayout(plotOutput("Expression"),plotOutput("Pseudotime")),
			 fixedRow(column(3, 
					 downloadButton('downloadPlot3')),
				  column(3, offset=3,
					 downloadButton('downloadPlot4'))))
	                )
	  )
)
)
)

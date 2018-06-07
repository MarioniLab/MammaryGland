library(shiny)
library(viridis)

shinyServer(

function(input, output) {

output$tSNE <- renderPlot({
	tsnPlot <- ggplot(pD, aes_string(x="tSNE1", y="tSNE2", color=input$colorBy)) +
	    geom_point(size=1) +
	    scale_color_manual(values=levels(pD$Colors)) +
	    theme_void()
	if (input$colorBy=="Condition") {
	    tsnPlot <- tsnPlot %+% scale_color_brewer(palette="Paired")
	}
	tsnPlot
    })

output$tSNE2 <- renderPlot({
	gene <- input$Gene
	gn <- filter(fD, symbol %in% gene) %>% .$id
	expr <- t(m[gn,])
	pD[,gn] <- as.vector(expr)
	pD <- arrange(pD, (pD[,gn]))
	tsnPlot <- ggplot(pD, aes_string(x="tSNE1", y="tSNE2", color=gn)) +
	    geom_point(size=1.5) +
	    scale_color_viridis() +
	    theme_void()
	tsnPlot
    })

output$distPlot1 <- renderPlot({
	gene <- input$Gene
	     plotGeneDist(m, pD, fD, gene, colorBy="Cluster")
       })

output$downloadPlot2 <- downloadHandler(
    filename <- function() { paste(input$Gene, '.pdf', sep='') },
    content <- function(file) {
	p <- plotGeneDist(m, pD, fD, input$Gene, colorBy="Cluster")
	pdf(file)
	print(p)
	dev.off()
    }
)

output$downloadPlot1 <- downloadHandler(
	filename <- function() { paste(input$Gene, '.pdf', sep='')},
	content <- function(file) {
	gene <- input$Gene
	gn <- filter(fD, symbol %in% gene) %>% .$id
	expr <- t(m[gn,])
	pD[,gn] <- as.vector(expr)
	pD <- arrange(pD, (pD[,gn]))
	tsnPlot <- ggplot(pD, aes_string(x="tSNE1", y="tSNE2", color=gn)) +
	    geom_point(size=1.5) +
	    scale_color_viridis() +
	    theme_void()
	pdf(file,width=10,height=7)
	print(tsnPlot)
	dev.off()
	}
    )

output$DiffMap <- renderPlot({
    pD.sub <- pD[!is.na(pD$DC1),]
    cols <- levels(pD.sub$Colors)
    p.clust <- ggplot(pD.sub, aes(x=DC1,y=DC2, color=Cluster)) +
	geom_point(size=2, pch=20) +
	scale_color_manual(values=cols) +
	xlab("Component 1") +
	ylab("Component 2") +
	guides(color=FALSE) +
	theme_void()
    p.clust
}
)

output$Pseudotime <- renderPlot({
    plotPseudotimeTrend(input)
}
)

output$Expression <- renderPlot({
    plotTriangleExpression(input)
}
)

output$downloadPlot3 <- downloadHandler(
	filename <- function() { paste(input$Gene, '.pdf', sep='')},
	content <- function(file) {
	    p <- plotTriangleExpression(input)
	    pdf(file)
	    print(p)
	    dev.off()
	}
	)

output$downloadPlot4 <- downloadHandler(
	filename <- function() { paste(input$Gene, '.pdf', sep='')},
	content <- function(file) {
	    p <- plotPseudotimeTrend(input)
	    pdf(file)
	    print(p)
	    dev.off()
	}
	)
}
)

# load shiny library
library(shiny)

# set up back end
shinyServer(function(input, output) {
    # load the data
    amlData <- read.delim("data/shinyExampleData.tsv")
    
    # construct a plot to show the data
    library(ggplot2)
    output$scatterPlot <- renderPlot({
        p1 <- ggplot(amlData, aes_string(x=input$x_axis, y=input$y_axis, colour=input$point_color)) + geom_point()
        p1 <- p1 + xlab("Variant Allele Fraction") + ylab("Variant Allele Fraction")
        p1
    })
})
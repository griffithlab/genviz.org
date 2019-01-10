# load shiny library
library(shiny)

# set up back end
shinyServer(function(input, output) {
  # load the data
  amlData <- read.delim("data/shinyExampleData.tsv")
  
  # construct a plot to show the data
  library(ggplot2)
  output$scatterPlot <- renderPlot({
    p1 <- ggplot(amlData, aes(x=Skin_d42_I_vaf, y=MC_d0_clot_A_vaf)) + geom_point()
    p1
  })
  
  output$text <- renderText({"This is some text"})
})
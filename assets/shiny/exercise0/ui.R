
# load shiny library
library(shiny)

# set up front end
shinyUI(fluidPage(
  mainPanel(plotOutput("scatterPlot"),
            textOutput("text"))
))
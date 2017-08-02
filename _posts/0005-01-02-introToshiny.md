---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to shiny
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-01-02
---

Interactive graphics is an emerging area particularly within R, there are many libraries available to make these sorts of visualizations however most of these libraries are still nascent. In this sub-module we will give a brief overview of [shiny](https://shiny.rstudio.com/), a web application framework within R for building interactive web pages. Using shiny we will build a simple application to display are data using reactive data sets and ggplot.

### Install shiny

The [shiny](https://shiny.rstudio.com/) package is available on cran and is fairly easy to install using [install.packages()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/install.packages). Go ahead and install and load the package. The package has 11 example apps built in which can be viewed using the [runExample()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/runExample) function, we will be building our own app from scratch, but feel free to try out a few of these examples to get a feel for what shiny can do.
```R
# install and load shiny
install.packages("shiny")
library(shiny)

# list the built in shiny app examples
runExample()

# run one of these examples in Rstudio
runExample("06_tabsets")
```

What shiny is actually doing here is converting the [R](https://www.r-project.org/) code to [html pages](https://en.wikipedia.org/wiki/HTML) and serving those on a random port using the ip address 127.0.0.1 which is [localhost](https://en.wikipedia.org/wiki/Localhost) on most computers. In simplified terms these [html pages](https://en.wikipedia.org/wiki/HTML) are simply being hosted by your own computer. If you are in Rstudio your web application should have been opened automatically, however you can also view these with any modern web browser by going to the web address listed after calling [runExample()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/runExample). It should look something like this: http://127.0.0.1:4379.

### Structure of a shiny app

You may have noticed some R code split into two files named server.R and ui.R when running the example above. The basic code to run any shiny app is split into two parts the server and user interface. the server script is the [back end](https://en.wikipedia.org/wiki/Front_and_back_ends) of our shiny web app and contains the instructions to build the app. The user interface script is the [front end](https://en.wikipedia.org/wiki/Front_and_back_ends) and is essentially what a user views and interacts with. Both of these files should be in the same directory for the app to work properly. Go ahead and make a folder for our shiny app called "testApp" and put these two scripts in there. This is the bare minimum for a shiny app and will generate an empty web application. Make sure that your current working directory in R is set to the top level of "testApp", you can use [getwd()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/getwd) and [setwd()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/getwd) to print and set this respectively.
* ui.R

```R
# load shiny library
library(shiny)

# set up front end
shinyUI(fluidPage(
))
```
* server.R

```R
# load shiny library
library(shiny)

# set up back end
shinyServer(function(input, output) {
})
```

### loading data into shiny
Now that we've got a basic frame work up let's go ahead and load some data and answer a few questions. The data we will use is supplemental table 6 from the paper ["Comprehensive genomic analysis reveals FLT3 activation and a therapeutic strategy for a patient with relapsed adult B-lymphoblastic leukemia."](https://www.ncbi.nlm.nih.gov/pubmed/27181063). The data contains capture variant sequencing information from an adult AML patient from 11 samples of various cell populations and timepoints. You can download the table [here](http://genomedata.org/gen-viz-workshop/shinyExampleData.tsv). We can load this data into shiny as you would any other data in R just be sure to do this in the server.R script. For simplicity make a "data" directory in your app and place the data file there. Then add this to your server.R script to make the data available within the back-end of shiny.
* server.R

```R
# load shiny library
library(shiny)

# set up back end
shinyServer(function(input, output) {
    # load the data
    amlData <- read.delim("data/shinyExampleData.tsv")
})
```

### sending output to the front end

* ui.R

```R
# load shiny library
library(shiny)

# set up front end
shinyUI(fluidPage(
    mainPanel(plotOutput("scatterPlot"))
))
```
* server.R

```R
# load shiny library
library(shiny)

# set up back end
shinyServer(function(input, output) {
    # load the data
    amlData <- read.delim("data/shinyExampleData.tsv")

    # construct a plot to show the data
    library(ggplot2)
    output$plot <- renderPlot({
        p1 <- ggplot(amlData, aes(x=Skin_d42_I_vaf, y=MC_d0_clot_A_vaf)) + geom_point()
        return(p1)
    })
})
```

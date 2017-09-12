---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to shiny
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-04-01
---

Interactive graphics is an emerging area within R. There are many libraries available to make interactive visualizations, however most of these libraries are still quite new. In this sub-module we will give a brief overview of [shiny](https://shiny.rstudio.com/), a web application framework within R for building interactive web pages. Using shiny we will build a simple application to display our data using reactive data sets and ggplot.

### Install shiny

The [shiny](https://shiny.rstudio.com/) package is available on cran and is fairly easy to install using [install.packages()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/install.packages). Go ahead and install and load the package. The package comes with 11 example apps that can be viewed using the [runExample()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/runExample) function, we will be building our own app from scratch, but feel free to try out a few of these examples to get a feel for what shiny can do. Shiny also provides a nice [gallery of example applications](https://shiny.rstudio.com/gallery/) and even a genomics example plotting [cancer genomics data in a circos-style application](https://shiny.rstudio.com/gallery/genome-browser.html).
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

After checkiny it out, use the escape key to stop the shiny app.

### Structure of a shiny app

The basic code to run any shiny app is split into two parts the server (e.g., server.R) and user interface (e.g., ui.R). The server script is the [back end](https://en.wikipedia.org/wiki/Front_and_back_ends) of our shiny web app and contains the instructions to build the app. The user interface script is the [front end](https://en.wikipedia.org/wiki/Front_and_back_ends) and is essentially what a user views and interacts with. Both of these files should be in the same directory for the app to work properly. 

Go ahead and make a folder for our shiny app called "testApp".

Next create the following two scripts there: `ui.R` and `server.R`. This is the bare minimum for a shiny app and will generate an empty web application.

* Put the following in a file: `ui.R`
```R
# load shiny library
library(shiny)

# set up front end
shinyUI(fluidPage(
))
```

* Put the following in a file: `server.R`
```R
# load shiny library
library(shiny)

# set up back end
shinyServer(function(input, output) {
})
```

To view/test your app simply type the `runApp(port=7777)` command in your R/Rstudio terminal. For convenience in this tutorial, we have selected a specific port instead of letting shiny choose one randomly. 

Make sure that your current working directory in R is set to the top level of "testApp" where you put server.R and ui.R. 

You can use [getwd()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/getwd) and [setwd()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/getwd) to print and set this respectively. 
test
So far, all you should see is an empty web page at [http://127.0.0.1:7777](http://127.0.0.1:7777).

Example:
```R
getwd()
setwd("/Users/mgriffit/Desktop/testApp)
getwd()
runApp(port=7777)
```

### Loading data into the shiny back end (server)
Now that we've got a basic frame work up let's go ahead and load some data and answer a few questions. The data we will use is supplemental table 6 from the paper ["Comprehensive genomic analysis reveals FLT3 activation and a therapeutic strategy for a patient with relapsed adult B-lymphoblastic leukemia."](https://www.ncbi.nlm.nih.gov/pubmed/27181063). The data contains variant allele frequency (VAF) values from a targeted capture sequencing study of an adult AML patient with 11 samples of various cell populations and timepoints. You can download the table [here](http://genomedata.org/gen-viz-workshop/intro_to_shiny/shinyExampleData.tsv). For simplicity, make a "data" directory in your app and place the data file there. We can load this data into shiny as you would any other data in R. Just be sure to do this in the server.R script and place the code within the unamed function. Add the following to your server.R script to make the data available within the shiny server.
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

### Sending output to the shiny front end (UI)
Now that we have data let's make a quick plot showing the distribution of VAF for the normal skin sample (Skin_d42_I_vaf) in comparison to the initial tumor marrow core sample (MC_d0_clot_A_vaf) and send it to the app's user interface. We'll need to first create the plot on the back end (i.e. server.R). We can use any graphics library for this, but here we use [ggplot2](http://ggplot2.tidyverse.org/reference/). In order to be compatible with the shiny UI we call a *Render* function, in this case [renderPlot()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/renderPlot) which takes an expression (i.e. set of instructions) and produces a plot. The curly braces in [renderPlot()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/renderPlot) just contain the expression used to create the plot and are useful if the expression takes up more than one line. The [renderPlot()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/renderPlot) will do some minimal pre-processing of the object returned in the expression and store it to the list-like "output" object. Notice that in the UI.R file we have added a [mainPanel()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/mainPanel) which, as it sounds, is instructing the app to create a main panel on the user interface. Now that we have somewhere to display our plot we can link what was created on the back end to the front end. This is done with the *Output* family of functions, in this case our output is a plot generated by [renderPlot()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/renderPlot) and is stored in the list like output object as "scatterplot" created in the server.R file. We use [plotOutput()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/plotOutput) to provide this link to the front end and give the output ID, which is just the name of the object stored in the output-like list. Note that when providing this link the type of object created with a *Render* function must correspond to the *Output* function, in this example we use [renderPlot()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/renderPlot) and [plotOutput()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/plotOutput) but other functions exist for other data types such as [renderText()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/renderText) and [textOuput()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/textOutput).

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
    output$scatterPlot <- renderPlot({
        p1 <- ggplot(amlData, aes(x=Skin_d42_I_vaf, y=MC_d0_clot_A_vaf)) + geom_point()
        p1
    })
})
```

Once again, to view/test your app simply type the `runApp(port=7777)` command in your R/Rstudio terminal and go to [http://127.0.0.1:7777](http://127.0.0.1:7777). This should happen automatically from Rstudio. If your previous app is still running you may need to stop and restart it and/or refresh your browser. You should now see a ggplot graphic in your browser (see below). But, so far, nothing is interactive about this plot. We will allow some basic user input and interactivity in the next section.

{% include figure.html image="/assets/shiny/shiny_simple_output.png" width="1000" %}

### Sending input from the front end
Now that we know how to link output from the back end to the front end let's do the opposite and link user input from the front end to the back end. Essentially this is giving the user control to manipulate user interface objects. Specifically let's allow the user to choose which sample Variant Allele Fraction (VAF) columns in the data set to plot on the x and y axis of our scatter plot. Let's start with the ui.R file. Below, we have added the [sidebarLayout()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/sidebarLayout) schema which will create a layout with a side bar and a main panel. Within this layout we define a [sidebarPanel()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/sidebarPanel) and a [mainPanel()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/mainPanel). Within the [sidebarPanel()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/sidebarPanel) we define two drop down selectors with [selectInput()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/selectInput). Importantly, within these functions we assign an *inputId* which is what will be passed to the back end. On the back end side (server.R) we've already talked about output within the unnamed function, a second argument exists as well called "input". This is the argument used to communicate from the front end to the back end and in our case it holds the information passed from each [selectInput()]([selectInput()](https://www.rdocumentation.org/packages/shiny/versions/1.0.3/topics/selectInput)) call with the id's "x_axis" and "y_axis". To make our plot reactively change based on this input we simply call up this information within the ggplot call. You might have noticed that we are using [aes_string()](http://ggplot2.tidyverse.org/reference/aes_.html) instead of [aes()](http://ggplot2.tidyverse.org/reference/aes_.html). This is only necessary because "input$x_axis" and "input$y_axis" are passed as strings and as such we need to let ggplot know this so the non-standard evalutation typically used with [aes()](http://ggplot2.tidyverse.org/reference/aes_.html) is not performed.

* ui.R

```R
#load shiny library
library(shiny)

# define the vaf column names
axis_options <- c("Skin_d42_I_vaf", "MC_d0_clot_A_vaf", "MC_d0_slide_A_vaf", "BM_d42_I_vaf",
                  "M_d1893_A_vaf", "M_d3068_A_vaf", "SB_d3072_A_rna_vaf", "SB_d3072_A_vaf",
                  "BM_d3072_A_vaf", "SL_d3072_I_vaf", "MC_d3107_A_vaf", "BM_d3137_I_vaf",
                  "M_d3219_I_vaf", "BM_d4024_I_vaf")

# set up front end
shinyUI(fluidPage(

  # set up the UI layout with a side and main panel
  sidebarLayout(

    # set the side panel to allow for user input
    sidebarPanel(
      selectInput(inputId="x_axis", label="x axis", choices=axis_options, selected="Skin_d42_I_vaf"),
      selectInput(inputId="y_axis", label="y axis", choices=axis_options, selected="MC_d0_clot_A_vaf")
    ),

    # set the plot panel
    mainPanel(
      plotOutput("scatterPlot")
    )
  )
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
  output$scatterPlot <- renderPlot({
    p1 <- ggplot(amlData, aes_string(x=input$x_axis, y=input$y_axis)) + geom_point()
    p1 <- p1 + xlab("Variant Allele Fraction") + ylab("Variant Allele Fraction")
    p1
  })
})
```

Once again, to view/test your app simply type the `runApp(port=7777)` command in your R/Rstudio terminal and go to [http://127.0.0.1:7777](http://127.0.0.1:7777). This should happen automatically from Rstudio. If your previous app is still running you may need to stop and restart it and/or simply refresh your browser. You should now see a ggplot scatterplot graphic in your browser (see below) as before. But, now you should also see user-activated drop-down menus that allow you to select which data to plot and visualize. You have created your first interative shiny application!

{% include figure.html image="/assets/shiny/shiny_interactive_output.png" width="1000" %}

### Exercises
We have given a very quick overview of [shiny](https://shiny.rstudio.com/), and have really only scraped the surface of what [shiny](https://shiny.rstudio.com/) can be used for. Using the knowledge we have already learned however let's try modifying our existing shiny app.

Right now the plot looks fairly bland. Try adding the ability for the user to enter a column name as text to color points by. For example, try coloring by the column names "Class" or "Clonal.Assignment". Use your existing ui.R and server.R files as a starting point. If successful, you should be able to restart/refresh your shiny app and see something like the following:

{% include figure.html image="/assets/shiny/shiny_exercise_output.png" width="1000" %}

{% include question.html question="Get a hint!" answer='You will want to use <a href="https://www.rdocumentation.org/packages/shinybootstrap2/versions/0.2.1/topics/textInput">textInput()</a> within the ui.R file for this and then link the input to the ggplot call.'%}

{% include question.html question="Solution" answer='These files contain the correct answer: <a href="http://genviz.org/assets/shiny/exercise1/ui.R">ui.R</a>, <a href="http://genviz.org/assets/shiny/exercise1/server.R">server.R</a>'%}

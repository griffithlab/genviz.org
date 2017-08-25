---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to transition/transverion plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-01
---

Transitions and transversions describe the mutation of a nucleotide from one base to another. Specifically transitions describe the mutation of nucleotides between purines or between pyrimidines and transversions from purines to pyrimidines or pyrimidines to purines. Often such events can give clues as to the origin of mutations within a sample. For example transverions are associated with ionizing radiation.

{% include figure.html image="/assets/GenVisR/Transitions_and_transversions.svg" width="350" link="https://commons.wikimedia.org/wiki/File:Transitions_and_transversions.svg" title="Transitions and Transversions" author="Krishnavedala" license="CC BY-SA 4.0" license_link="https://creativecommons.org/licenses/by-sa/4.0/deed.en" position="center" %}

### Loading Data
Figure 1 of the manuscript ["Recurrent somatic mutations affecting B-cell receptor signaling pathway genes in follicular lymphoma"](http://www.bloodjournal.org/content/129/4/473/tab-figures-only?sso-checked=true) plots the rate of each type of transition and transversion observed in a cohort of samples. In this section we will use the [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) function [TvTi](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) to create a version of this figure, shown below. The underlying data for this figure can be found in supplemental table S5 and is the same data used in the ["Introduction to ggplot2"](http://localhost:4000/module%202/0002/03/01/introToggplot2/) section. You can download the data  [here](http://www.genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv) if you haven't already.

{% include figure.html image="/assets/GenVisR/Figure1_Discovery_all_mutrate_tvti_v2.png" width="650" link="https://commons.wikimedia.org/wiki/File:Transitions_and_transversions.svg" title="This research was originally published in Blood. Krysiak et al. Blood 2017 129:473-483" author="Krysiak et al." license="&copy; the American Society of Hematology" license_link="http://www.bloodjournal.org/content/129/4/473/tab-figures-only?sso-checked=true" %}

### Data Preparation

The [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) function [TvTi()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) will calculate the rate of individual transitions and transversions within a cohort and plot this information as a stacked barchart. As with the [waterfall()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/waterfall) function [TvTi()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) can accept multiple file types as specified by the parameter `fileType`. Here we will use the most flexible, an "MGI" file type which expects a data frame with column names "sample", "reference", and "variant". Further we will re-assign the sample column to use the same sample names as in the manuscript figure. In order to match the manuscript figure we will also limit the data to just the discovery cohort and re-factor the data frame with the function [factor()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/factor) so samples without any data are not plotted.

```R
# load the data into R
mutationData <- read.delim("ggplot2ExampleData.tsv")

# change the "Simple_name" column to "sample"
colnames(mutationData)[colnames(mutationData) %in% "sample"] <- "sample_1"
colnames(mutationData)[colnames(mutationData) %in% "Simple_name"] <- "sample"

# subset the data to just the discovery cohort
mutationData <- mutationData[mutationData$dataset == "discovery",]
mutationData$sample <- factor(mutationData$sample, levels=unique(mutationData$sample))

# run TvTi
TvTi(mutationData, fileType="MGI")
```

{% include figure.html image="/assets/GenVisR/transition_transversion_v2.png" width="750" %}

### Plotting Frequencies with Proportions

You might notice from the previous plot that two samples look a bit off. Specifically the samples "LYM002-Naive" and "LYM005-Naive" are missing transitions and transversions. By default [TvTi()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) plots only the proportion of transitions and transversions which can be deceiving if the frequency of mutation events is not high enough to get an accurate calculation of the relative proportion. We can use the parameter `type` which expects a string specifying either "proportion" or "frequency", to instead plot the frequency of each mutation event. Let's use this parameter two create two plots, one of frequency the other of proportion and combine them. To do this we will need to output our plots as grobs using `out="grob"`, and use the function [arrangeGrob()](https://www.rdocumentation.org/packages/gridExtra/versions/2.2.1/topics/arrangeGrob) from the [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) package to combine the two plots into one. Finally we will need to call [grid.draw()](https://www.rdocumentation.org/packages/grid/versions/3.4.1/topics/grid.draw) on the grob to output it to the current graphical device. We can see from the resulting plot that quite a few samples have low mutation frequencies and we should interpret these samples with caution.

```R
# install and load gridExtra
install.packages("gridExtra")
library(gridExtra)
library(grid)

# make a frequency and proportion plot
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency")
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion")

# combine the two plots
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1)

# draw the plot
grid.draw(finalGrob)
```

{% include figure.html image="/assets/GenVisR/transition_transversion_v3.png" width="750" %}


### Changing the order of samples

By default the order of the samples on the x-axis uses the levels of the sample column in the input data frame to determine an order, if that column is not of type factor it is coerced to one using the existing order of samples to define the levels within that factor. This behavior can be altered using parameter `sort` which will take a character vector specifying either "sample" or "tvti". Specifying "sample" will order the x-axis alpha-numerically based on the sample names, conversely "tvti" will attempt to order the x-axis samples by decreasing proportions of each transition/transversion type. In the manuscript figure we can see that a custom sort order is being used. Let's take advantage of the default ordering behavior in [TvTi()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) and change the levels of the sample column to match the order in the manuscript figure.

```R
# determine the exisitng order of samples
levels(mutationData$sample)

# create a custom sample order to match the manuscript
sampleOrder <- c("LYM145-Naive", "LYM023-Naive", "LYM045-Naive", "LYM238-Sorted", "LYM238-Naive", "LYM058-Naive", "LYM125-Naive", "LYM100-Naive", "LYM088-Naive", "LYM167-Naive", "LYM153-Naive", "LYM002-Naive", "LYM005-Naive", "LYM120-Naive", "LYM120-Treated", "LYM215-Sorted", "LYM215-Treated", "LYM139-Treated", "LYM177-Treated", "LYM202-Sorted", "LYM202-Treated", "LYM062-Treated", "LYM057-tNHL", "LYM001-tNHL", "LYM175-tNHL", "LYM193-tNHL", "LYM237-tNHL", "LYM013-tNHL")

# make sure all samples have been specified
all(mutationData$sample %in% sampleOrder)

# alter the sample levels of the input data
mutationData$sample <- factor(mutationData$sample, levels=sampleOrder)

# recreate the plot
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency")
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion")
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1)
grid.draw(finalGrob)
```

{% include figure.html image="/assets/GenVisR/transition_transversion_v4.png" width="750" %}

### adding clinical data

As with many other [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) functions we can add a heatmap of clinical data via the parameter `clinData` which expects a data frame with column names "sample", "variable", and "value". A subset of the clinical data is available to download [here](http://genomedata.org/gen-viz-workshop/GenVisR/FL_ClinicalData.tsv), go ahead and load it into R. The clinical data contains information for all samples in the experiment however we are only interested in the discovery cohort, the next step then is to subset our clinical data to only those samples in our main plot. From there we can use [melt()](https://www.rdocumentation.org/packages/reshape2/versions/1.4.2/topics/melt) to coerce the data frame into the required "long" format and add the clinical data to the proportions plot with the `clinData` parameter. In order to more closely match the manuscript figure we will also define a custom color pallette and set the number of columns in the legend with the parameters `clinVarCol` and `clinLegCol` respectively.

```R
# read in the clinical data
clinicalData <- read.delim("FLClinicalData.tsv")

# subset just the discovery cohort
clinicalData <- clinicalData[clinicalData$Simple_name %in% mutationData$sample,]
all(sort(unique(as.character(clinicalData$Simple_name))) == sort(unique(as.character(mutationData$sample))))

# convert to long format
library(reshape2)
clinicalData <- melt(data=clinicalData, id.vars="Simple_name")

# define a few clinical parameter inputs
clin_colors <- c('0'="lightblue",'1'='dodgerblue2','2'='blue1','3'='darkblue','4'='black','FL'='green3','t-NHL'='darkgreen','Treated'='darkred','TreatmentNaive'='red','NA'='lightgrey')

# add the clinical data to the plots
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency")
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion", clinData=clinicalData, clinLegCol = 2, clinVarCol = clin_colors)
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1, heights=c(2, 5))
grid.draw(finalGrob)
```

{% include figure.html image="/assets/GenVisR/transition_transversion_v5.png" width="750" %}

### re-aligning plots

Up til this point we have glossed over what [arrangeGrob()](https://www.rdocumentation.org/packages/gridExtra/versions/2.2.1/topics/arrangeGrob) is actually doing however this is no longer a luxury we can afford. In the previous figure we can see that adding the clinical data to our plot has caused our proportion and frequency plots to go out of alignment. This has occurred because the y-axis text in the clinical plot takes up more space than in the frequency plot. Because we added clinical data with [GenVisR]() to the proportion plot those graphics are aligned with each other however [TvTi()]() has no idea the frequency plot exists, it's added to our graphic after the fact with [arrangeGrob()](). In order to fix this issue we will need a basic understanding of "grobs", "TableGrobs", and "viewports". First off a [grob]() is just short for "grid graphical object" from the low-level graphics package [grid](); Think of it as a set of instructions for create a graphical object (i.e. a plot). The graphics library underneath all of [ggplot2's]() graphical elements are really composed of [grob's]() because ggplot2 uses [grid]() underneath. A  TableGrob is a class from the [gtable]() package and provides an easier way to view and manipulate groups of [grobs](), it is actually the intermediary between [ggplot2]() and [grid](). A "viewport" is a graphics region for which a [grob]() is assigned. We have already been telling [GenVisR]() to output [grobs]() instead of drawing a graphic because we have been using the [arrangeGrob()]() function in [gridExtra]() to arrange and display extra viewports in order to create a final plot. Let's briefly look at what this actually means, we can show the layout of viewports with the function [gtable_show_layout](). In the figure below we have overlayed this layout on to our plot and have shown the gtable to the right. In the table grob we can see that at the top level our [grobs]() are compsed of 2 vertical elements and 1 horizontal elements. The column z corresponds to the order of plot layers, and cells corresponds to the viewport coordinates. So we can see that the [grob]() listed in row 1 is plotted first and is located on the viewport spanning vertical elemnents 1-1 and horizontal elements 1-1 and thus on our layout is in the viewport (1, 1).

```R
# load the gtable library
library(gtable)

# visualize the viewports of our plot
gtable_show_layout(finalGrob)
```

{% include figure.html image="/assets/GenVisR/transition_transversion_v6.png" width="750" %}

Pretty simple right? But not so fast the TableGrob in `finalGrob` is a little over simplified on the surface. It is true that our plot is composed of two final [grobs]() however each [grob]() is made up of lists of [grobs]() in the GrobTable all the way down to the individual base elements which make up the plot. In the figure below we illustrate the second layer of these lists. To re-align the plots we'll need to peel back the layers until we can acess the original [grobs]() from [ggplot2]. We will then need to obtain the widths of each of these elements, find the max width with [unit.pmax](), and re-assign the width of each as these [grobs]() using the max width.

{% include figure.html image="/assets/GenVisR/transition_transversion_v7.png" width="750" %}

```R
# find the layer of grobs widths we care about
proportionPlotWidth <- finalGrob$grobs[[2]]$grobs[[1]]$widths
clinicalPlotWidth <- finalGrob$grobs[[2]]$grobs[[2]]$widths
transitionPlotWidth <- finalGrob$grobs[[1]]$grobs[[1]]$widths

# find the max width of each of these
maxWidth <- unit.pmax(proportionPlotWidth, clinicalPlotWidth, transitionPlotWidth)

# re-assign the max widths
finalGrob$grobs[[2]]$grobs[[1]]$widths <- as.list(maxWidth)
finalGrob$grobs[[2]]$grobs[[2]]$widths <- as.list(maxWidth)
finalGrob$grobs[[1]]$grobs[[1]]$widths <- as.list(maxWidth)

# plot the plot
grid.draw(finalGrob)
```

### adding additional layers

Our plot is getting getting close to looking like the one in the manuscript, there are a few final touches we should add however. Specifically let's add in some vertical lines separating samples by "treatment-naive", "treated", and "t-NHL". Let's also remove the redundant x-axis text, we can do all this by adding additional [ggplot2]() layers to the various plots with the parameter `layers` and `clinLayer`.

```R
# load ggplot2
library(ggplot2)

# create a layer for vertical lines
layer1 <- geom_vline(xintercept=c(14.5, 22.5))

# create a layer to remove redundant x-axis labels
layer2 <- theme(axis.text.x=element_blank(), axis.title.x=element_blank())

# create the plot
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency", layers=list(layer1, layer2))
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion", clinData=clinicalData, clinLegCol = 2, clinVarCol = clin_colors, layers=layer1, clinLayer=layer1)
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1, heights=c(2, 5))
grid.draw(finalGrob)
```

{% include figure.html image="/assets/GenVisR/transition_transversion_v8.png" width="750" %}

#### Exercises

The plot above has been finnished, however we re-created all the plots without taking the time to re-align them. Try and fix that now, if you get stuck go over the aligning grobs section again.

{% include question.html question="Re-align the plot above." answer='Go through the aligning grobs section if you get stuck'%}

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

As with many other [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) functions we can add a heatmap of clinical data via the parameter `clinData` which expects a data frame with column names "sample", "variable", and "value".

### adding additional layers

###

```R
library(ggplot2)

layer1 <- geom_vline(xintercept=5)

# recreate the plot
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency", layers=layer1)
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion", layers=layer1)
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1)
grid.draw(finalGrob)
```

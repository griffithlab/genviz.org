---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to transition/transverion plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-01
---

Transitions and transversions describe the mutation of a nucleotide from one base to another. Specifically transitions describe the mutation of nucleotides between purines or between pyrimidines and transversions from purines to pyrimidines or pyrimidines to purines. Often such events can give clues as to the origin of mutations within a sample. For example transverions are associate with ionizing radiation.

{% include figure.html image="/assets/GenVisR/Transitions_and_transversions.svg" width="350" link="https://commons.wikimedia.org/wiki/File:Transitions_and_transversions.svg" title="Transitions and Transversions" author="Krishnavedala" license="CC BY-SA 4.0" license_link="https://creativecommons.org/licenses/by-sa/4.0/deed.en" position="center" %}

### Loading Data
Figure 1 of the manuscript ["Recurrent somatic mutations affecting B-cell receptor signaling pathway genes in follicular lymphoma"](http://www.bloodjournal.org/content/129/4/473/tab-figures-only?sso-checked=true) plots the rate of each type of transition and transversion. In this section we will use the [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) function [TvTi](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) to create a version of this figure, shown below. The underlying data for this figure can be found in supplemental table S5 and is the same data used in the ["Introduction to ggplot2"](http://localhost:4000/module%202/0002/03/01/introToggplot2/) section. You can download the data  [here](http://www.genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv) if you haven't already.

{% include figure.html image="/assets/GenVisR/Figure1_Discovery_all_mutrate_tvti_v2.png" width="650" link="https://commons.wikimedia.org/wiki/File:Transitions_and_transversions.svg" title="This research was originally published in Blood. Krysiak et al. Blood 2017 129:473-483" author="Krysiak et al." license="&copy; the American Society of Hematology" license_link="http://www.bloodjournal.org/content/129/4/473/tab-figures-only?sso-checked=true" %}

### Data Preparation

The [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) function [TvTi()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) will calculate the rate of individual transitions and transversions within a cohort and plot this information as a stacked barchart. As with the [waterfall()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/waterfall) function [TvTi()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) can accept multiple file types as specified by the parameter `fileType`. Here we will use the most flexible, an "MGI" file type which expects a data frame with column names "sample", "reference", and "variant". Further we will re-assign the sample column to use the same sample names as the manuscript figure. In order to match the manuscript figure we will also limit the data to just the discovery cohort and re-factor the data frame with the function [factor()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/factor) so samples without any data are not plotted.

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

You might notice from the previous plot that two samples look a bit off. Specifically the samples "LYM002-Naive" and "LYM005-Naive" are missing transitions and transversions. By default [TvTi()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/TvTi) plots only the proportion of transitions and transversions which can be deceiving if the frequency of mutation events is not high enough to get an accurate calculation of the relative proportion. We can use the parameter `type` which expects a string specifying either "proportion" or "frequency", to instead plot the frequency of each mutation event. Let's use this parameter two create two plots, one of frequency the other of proportion and combine them. To do this we will need to output our plots as grobs using `out="grob"` and use the function [arrangeGrob()](https://www.rdocumentation.org/packages/gridExtra/versions/2.2.1/topics/arrangeGrob) from the [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) package to combine the two plots into one. Finally we will need to call [grid.draw()](https://www.rdocumentation.org/packages/grid/versions/3.4.1/topics/grid.draw) on the grob to output it to the current graphical device. We can see from the resulting plot that quite a few samples have low mutation frequencies and we should interpret these samples with caution.

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

---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Advanced ggplot2
categories:
    - Module-02-R
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-02
---

We've gone over the basics of ggplot2 in the previous section, in this section we will go over some more advanced topics related to ggplot2 and its underlying concepts. We will explore how to modify core elements of a plot after it's created, how to plot separate plots on the same page, and how to make sure plots align to one another.

#### Aligning plots on the same page

Often when creating figures for publications the figure in the final manuscript is actually composed of multiple figures. You could of course achieve this using third-party programs such as illustrator and inkscape however that means having to manually redo the figure in those programs anytime there is an update or change. Instead of going through that hassel let's learn how to do it all programatically. To achieve this we will of course make our plots using ggplot, and use the gridExtra package to do the actual aligning. Let's go ahead and get started.

Our first step is to make some plots. We'll go ahead and use the same dataset as before which is available at [http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv](http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv) as a reminder. Let's go ahead and load this data in, if it isn't in R already.

```R
# install the ggplot2 library and load it
library(ggplot2)

# load Supplemental Table S5
# note that in the following example we are loading directly from a URL (instead of dowloading it to the instance first)
variantData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv")
```

Next let's start by making a bar chart for the genes which have over 20 mutations in all samples. The commands here should be review from previous sections.

```R
# install the plyr library
library(plyr)

# count the number of variants for each gene/mutation type and format for ggplot
geneCount <- count(variantData, vars=c("gene_name", "type"))
geneCount <- geneCount[geneCount$freq > 20,]

# to order samples for ggplot also get a total count overall
geneCountOverall <- aggregate(data=geneCount, freq ~ gene_name, sum)
geneOrder <- geneCountOverall[order(geneCountOverall$freq),]$gene_name
geneCount$gene_name <- factor(geneCount$gene_name, levels=geneOrder)

# make the barchart
p1 <- ggplot() + geom_bar(data=geneCount, aes(x=gene_name, y=freq, fill=type), stat="identity") + xlab("Gene") + ylab("Frequency") + scale_fill_manual("Mutation", values=c("#F97F51", "#55E6C1")) + theme_bw() + theme(plot.background = element_rect(color="red", size=2))

```

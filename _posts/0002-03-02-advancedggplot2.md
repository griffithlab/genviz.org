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

We can see that KMT2D seems interesting, it's got the most mutations and a fair number of deletions as well. To explore this a bit further lets go ahead and compare KMT2D to a few of the other highly conseverved genes from our first plot.

```R
geneCompare1 <- variantData[variantData$gene_name %in% c("KMT2D", "TNFRSF14"),]
geneCompare1 <- geneCompare1[,c("gene_name", "trv_type", "ucsc_cons")]
geneCompare1 <- geneCompare1[geneCompare1$trv_type == "missense",]
geneCompare1$ucsc_cons <- as.numeric(as.character(geneCompare1$ucsc_cons))
p2 <- ggplot() + geom_boxplot(data=geneCompare1, aes(x=gene_name, y=ucsc_cons, fill=gene_name)) + scale_fill_manual("Gene", values=c("#e84118", "#e1b12c")) + theme_bw() + xlab("Gene") + ylab("Conservation\nscore") + theme(plot.background = element_rect(color="dodgerblue", size=2))
```

```R
geneCompare2 <- variantData[variantData$gene_name %in% c("KMT2D", "BCL2"),]
geneCompare2 <- geneCompare2[,c("gene_name", "trv_type", "ucsc_cons")]
geneCompare2 <- geneCompare2[geneCompare2$trv_type == "missense",]
geneCompare2$ucsc_cons <- as.numeric(as.character(geneCompare2$ucsc_cons))
p3 <- ggplot() + geom_boxplot(data=geneCompare2, aes(x=gene_name, y=ucsc_cons, fill=gene_name)) + scale_fill_manual("Gene", values=c("#e84118", "#4cd137")) + theme_bw() + xlab("Gene") + ylab("Conservation\nscore") + theme(plot.background = element_rect(color="green", size=2))
```

```R
# convert to grobs
library(gridExtra)
p1_grob <- ggplotGrob(p1)
p2_grob <- ggplotGrob(p2)
p3_grob <- ggplotGrob(p3)

# make a layout
layout <- rbind(c(1, 1),
                c(2, 3))
grid.arrange(p1_grob, p2_grob, p3_grob, layout_matrix=layout)

```

```R
p4 <- ggplot() + geom_bar(data=geneCompare1, aes(x=gene_name, fill=gene_name)) + scale_fill_manual("Gene", values=c("#e84118", "#e1b12c")) + theme_bw() + theme(plot.background = element_rect(color="darkorange2", size=2)) + xlab("Frequency") + ylab("Gene")
p5 <- ggplot() + geom_bar(data=geneCompare2, aes(x=gene_name, fill=gene_name)) + scale_fill_manual("Gene", values=c("#e84118", "#4cd137")) + theme_bw() + theme(plot.background = element_rect(color="black", size=2)) + xlab("Frequency") + ylab("Gene")
```

```R
p4_grob <- ggplotGrob(p4)
p5_grob <- ggplotGrob(p5)

layout <- rbind(c(1, 1),
                c(2, 3),
                c(4, 5))
grid.arrange(p1_grob, p4_grob, p5_grob, p2_grob, p3_grob, layout_matrix=layout)
```

```R
# align plots
library(grid)
p4_grob_widths <- p4_grob$widths
p5_grob_widths <- p5_grob$widths
p2_grob_widths <- p2_grob$widths
p3_grob_widths <- p3_grob$widths

maxWidth <- unit.pmax(p4_grob_widths, p5_grob_widths, p2_grob_widths, p3_grob_widths)

p4_grob$widths <- maxWidth
p5_grob$widths <- maxWidth
p2_grob$widths <- maxWidth
p3_grob$widths <- maxWidth

layout <- rbind(c(1, 1),
                c(2, 3),
                c(4, 5))
grid.arrange(p1_grob, p4_grob, p5_grob, p2_grob, p3_grob, layout_matrix=layout)
```

```R
p2 <- p2 + theme(legend.position="none")
p3 <- p3 + theme(legend.position="none")

p2_grob <- ggplotGrob(p2)
p3_grob <- ggplotGrob(p3)

p2_grob <- gtable_add_cols(p2_grob, widths=unit(1, "null"), pos=8)
p2_grob <- gtable_add_cols(p2_grob, widths=unit(1, "null"), pos=8)

p3_grob <- gtable_add_cols(p3_grob, widths=unit(1, "null"), pos=8)
p3_grob <- gtable_add_cols(p3_grob, widths=unit(1, "null"), pos=8)

maxWidth <- unit.pmax(p4_grob_widths, p5_grob_widths, p2_grob_widths, p3_grob_widths)

p4_grob$widths <- maxWidth
p5_grob$widths <- maxWidth
p2_grob$widths <- maxWidth
p3_grob$widths <- maxWidth

layout <- rbind(c(1, 1),
                c(2, 3),
                c(4, 5))
finalGrob <- grid.arrange(p1_grob, p4_grob, p5_grob, p2_grob, p3_grob, layout_matrix=layout)
```

```R
finalGrob$grobs[[1]]$grobs[[7]]$children$axis$grobs[[2]]$children$GRID.text.6880$gp$col <- c("blue", "blue", "blue", "blue", "blue", "red", "red", "blue", "red")
```

---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to Loss of Heterozygosity Plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-07-01
---

We've gone through visualizations of point mutations and copy number changes using [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html), another type of genomic alteration that is often usefull to visualize are ["Loss of Heterozygosity"](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) events. In a diploid cell there are pairs of chromosomes each containing a single copy of the genome. These pairs each come from haploid gametes and are slightly different from each other leading to heterozygosity throughout much of the genome. Situations can arise however where this inherit heterozygosity of the genome is loss, commonly this is through deletions of a parental copy within a chromosome region also known as hemizygosity. Viewing deletions however will not give a complete picture of [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) as events can arise that will lead to copy-neutral LOH, for example if a parental copy was deleted but the then remaining copy underwent an amplification. In this section we will use the [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) function [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) created specifically for viewing loh events.

### How lohSpec works
The function [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) works by using a sliding window approach to calculate the difference in the variant allele fractions (VAF) between heterozygous variants in matched tumor normal samples. Essentially, a window of a specified size (defined by the parameter `window_size`) will obtain the absolute difference between the tumor VAF and the normal VAF, which is assumed to be .5, for each genomic position within the window's position. All of these points are then averaged to obtain a measure of [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) within the window. This window will then move forward by a length specified by the parameter `step`, and again cauclate this absolute mean tumor-normal VAF difference metric. This [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) value across all overlapping windows is averaged.

### Introduction to demonstration dataset
For this section we will be starting from a file of variants generated from the variant calling algorithm [varscan](http://varscan.sourceforge.net/). These variants originate from the breast cancer cell line [HCC1395](https://www.atcc.org/Products/All/CRL-2324.aspx) aligned to "hg19" and were lightly processed to limit to only events called as "Germline" or "LOH" by varscan. You can download the variant file from this [link](http://genomedata.org/gen-viz-workshop/GenVisR/HCC1395.varscan.tsv). Once you have the file we'll need to do some minor reformatting as described below to make it compatiable with [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec).

```R
# read in the varscan data
lohData <- read.delim("~/Desktop/HCC1395.varscan.tsv", header=FALSE)

# grab only those columns which are required and name them
lohData <- lohData[,c("V1", "V2", "V7", "V11")]
colnames(lohData) <- c("chromosome", "position", "n_vaf", "t_vaf")

# add a sample column
lohData$sample <- "HCC1395"

# convert the normal and tumor vaf columns to fractions
lohData$n_vaf <- as.numeric(gsub("%", "", lohData$n_vaf))/100
lohData$t_vaf <- as.numeric(gsub("%", "", lohData$t_vaf))/100

# limit to just the genome of interest
lohData <- lohData[grepl("^\\d|X|Y", lohData$chromosome),]
```

The function [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) requires a data frame of heterozygous variant calls with column names "chromosome", "position", "n_vaf", "t_vaf", and "sample" as it's primary input. Above we read in our germline variant calls from [varscan](http://varscan.sourceforge.net/) and subset the data to just the required columns renaming them with [colnames()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/row%2Bcolnames) and addin in a "sample" column. If we look at the the documentation we can see that [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) expects proportions for the VAF columns and not percentages. We coerce these columns into this format by using [gsub()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/grep) to remove the "%" symbol and dividing the resultant numeric value from the call [as.numeric()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/numeric) by 100. We also subset our data using [grepl()]() to remove anything that is not a canonical chromosome as the chromosomes in the primary data must match those specified in the `genome` parameter we will use later.

# Creating an initial plot
Like it's counterpart [cnSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/cnSpec), [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) needs genomic boundaries to be specified in addition to the primary data in order to ensure that the entire genome is plotted when creating the graphic. This can be done via the parameter `genome` which expects a character string specifying one of "hg19", "hg38", "mm9", "mm10", or "rn5". Alternatively a custom genome can be supplied with the `y` parameter. These parameters are both identical to their [cnSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/cnSpec) counterparts, please refer to the cnSpec [tutorial](http://genviz.org/module%203/0003/06/01/cnSpec_GenVisR/) for a more detailed explanation. Now that we have all of the required data let's go ahead and make a standard [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) plot.

```R
# Create an inital plot
lohSpec(lohData, genome="hg19")
```

{% include figure.html image="/assets/GenVisR/lohSpec_v1.png" width="950" %}

We have our first plot however you might have noticed a warning message along the lines of "Detected values with a variant allele fraction either above .6 or below .4 in the normal. Please ensure variants supplied are heterozygous in the normal!". Admittedly this warning message is accurate, let's take a minute to consider why biologically [lohSpec()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) is producing this warning message.

{% include question.html question="Get a hint!" answer='Refer to the "How lohSpec works" section above' %}

{% include question.html question="Answer" answer='A heterozygous call in a diploid genome should have a VAF of .5, any deviation from this is either noise due to insufficient coverage, or more problematically indicative that the site is not completely heterozygous. A reasonable explanation could be a focal amplification at that site on one allele. In any case it is inadvisable to use such sites when constructing this type of plot' %}

Now let's fix the issue and reproduce our plot, note that as expected the warning message has disappeared.

```R
# Obtain variants with a VAF greater than 0.4 and below 0.6.
lohData <- lohData[lohData$n_vaf > 0.4 & lohData$n_vaf < 0.6,]

# run lohSpec
lohSpec(x=lohData)
```

{% include figure.html image="/assets/GenVisR/lohSpec_v1.png" width="950" %}

Now that we have a plot let's go over what it is showing. First off what is a high [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) value. Remember the formula used to calculate [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) for a single variant is `|normal vaf - tumor vaf|` and that these values are averaged for a given window. Consider this example of [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) then `|.5 - 1 | = .5` and `|.5 - 0| = .5`. Conversely consider this example of non [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity), `| .5 - .5 | = 0`. So we can see that the closer we get to the value `.5` the more evidence there is that an [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) event exists. We have to take noise into account but we can see that the majority of this genome has some sort of [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity), probably not surprising as this data is derived from an imortalized cell line.

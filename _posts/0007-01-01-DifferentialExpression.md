---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Differential Expression with DEseq2
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0007-01-01
---

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Oftentimes, it will be used to define the differences between multiple biological conditions (e.g. post-treatment vs. untreated samples). There are many tools available to perform this type of analysis, in this course we will rely on the bioconductor package [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). We will then make various visualizations to help interpret our results.

### dataset
For this analysis we will use the RNAseq data from [E-GEOD-50760](https://www.ncbi.nlm.nih.gov/pubmed/25049118). This data consists of 54 samples from 18 individuals, each individual has a primary colorectal cancer sample, a metastatic liver sample, and a normal sample of the surrounding colonic epithilium. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts. We will use the output from HTseq as a starting point which can be downloaded [here](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts), sample information can be downloaded [here](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/ExperimentDesignFile.RnaSeq/experiment-design). A full description of the experimental design can be found at [array express](http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-50760/) and the [expression atlas](http://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Results).

### How does DEseq2 work
[DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a popular differential expression analysis package available through bioconductor. It's differential expression tests are based on a negative binomial generalized linear model. To get started we will first need to install the package and load the library so let's begin.
```R
# install the latest version of DEseq2
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
```
### Input data
Input data for [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) consists of un-normalized sequence read counts at either the gene or transcript level. No normalization of this data is needed because [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) internally corrects for factors, specifically library size using these raw counts. The tool [HTseq](http://htseq.readthedocs.io/en/release_0.9.0/) can be used to obtain this information and is what was used for our [example data](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts), let's go ahead and load this data and the [sample information](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/ExperimentDesignFile.RnaSeq/experiment-design) into R.

```R
# Read in the raw read counts
rawCounts <- read.delim("E-GEOD-50760-raw-counts.tsv")

# Read in the sample mappings
sampleData <- read.delim("E-GEOD-50760-experiment-design.tsv")
```

The next step is to create an object of class DESeqDataSet, this will store the readcounts and intermediate calculations needed for the differential expression analysis. The object will also store the design formula which is used to estimate dispersion and log2 fold changes used within the model. When specifying the formula it should take the form of a ~ followed by + signs separating factors. When using the default DEseq2 parameters the factor of interest (tissue type in this case) should be specified last and the control within that factor should be first when viewing the [levels()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/levels) for that variable. There are 4 methods to create this object depending on the format the input data is in. Because we already have our data loaded into R we will use [DESeqDataSetFromMatrix()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class).

```R
# convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID

# convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.clinical.information.", "Sample.Characteristic.individual.")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("tissueType", "individualID")
sampleData$individualID <- factor(sampleData$individualID)

# put the columns of the count data in the same order as rows names of the sample mapping
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# Check the levels of tissue type to make sure the control variable is first
sampleData$tissueType <- factor(sampleData$tissueType, levels=c("normal-looking surrounding colonic epithelium", "primary colorectal cancer", "metastatic colorectal cancer to the liver"))

# create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ individualID + tissueType)
```

This was quite a bit of code, let's go over whats going on here. The first thing we do is coerce the data frame containing the read counts into a format DESeq2 can accept. Specifically this must be a matrix with row names as genomic features (i.e. genes), and column names as samples. Next DESeq2 requires a data frame specifying the mapping of samples to variables, we load this in and clean it up some keeping only the variables we care about and making sure everything is a factor. For DEseq2 to work properly the column names of the count matrix must be in the same order as the row names of the sample mapping data, to ensure this we re-order the column names of the count data and run a check to ensure this has occurred correctly. To take advantage of the default settings of DEseq2 the control of the variable of interest, in our case the tissue type, should be the first element in the levels of that variable. Because we have more than 2 conditions for this variable we will not be taking advantage of the default settings however it's good to get into the practice of doing this so we do it here. We then create a DEseq2DataSet object with this information and supply a formula where we use the individual id as a blocking factor and tissue type as the comparison variable.

### Pre-filtering of data
While not strictly necessary it is good to do some preliminary filtering of the data before running the differential expression analysis. This will reduce the size of the DESeq2 object and speed up the the speed of the algorithm. Here we are performing relatively minor filtering requiring genes to have more than 5 read of support.
```R
# perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]
```

# set up multi-cores (optional)
The next two steps can take some time to perform, we can offset this somewhat by enabling multiple cores using [BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html). To take advantage of this you will need to install the [BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html) library and register the number of cores to use depending on your machine. Then when calling [DESeq()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeq) and [results()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results) add parallel=TRUE as a parameter to these function calls.
```R
# install and load the library
source("https://bioconductor.org/biocLite.R")
biocLite("BiocParallel")

# register the number of cores to use
register(MulticoreParam(4))
```

### Differential Expression Analysis
The next step is to run the function [DEseq()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeq) on our DESeq2 data set object. In this step the algorithm will perform the following steps:
1. estimation of size factors
2. estimation of dispersion
3. Negative Binomial GLM fitting and Wald statistic.

This step can take a few minutes to perform, for convenience a .RData object containing an R environment up to this step is available to download [here](http://genomedata.org/gen-viz-workshop/differentialExpression.RData). You can load this into your R environment with [load()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/load).
```R
# run pipeline for differential expression steps
deseq2Data <- DESeq(deseq2Data)

# load the R environment with this object (optional)
load(differentialExpression.RData)
```

### Extracting results
Finally we can extract the differential expression results with the [results()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results) function. When using this function we need to tell DESeq2 what comparison to make. This is only neccessary if the design formula is multi-factorial or as in our case the variable in the design formula has more than 2 levels. This is done with the contrast parameter which takes a character vector of three elements giving the name of the factor of interest, the numerator, and the denominator (i.e. control). Let's get output for normal vs primary expression results and view a summary of results.
```R
# Extract differential expression results
deseq2Results <- results(deseq2Data, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))

# view summary of results
summary(deseq2Results)
```
### MA-plot
An MA plot plots a log ratio (M) over an average (A), it is a way to visualize the differences between two groups. In general we would expect the expression of genes to remain consistent between conditions and so the MA plot should be similar to the shape of a trumpet with most points residing on the 0 intercept of the y axis. DEseq2 has a built in method for constructing an MA plot of our results however as this is a a visualization course let's go ahead and construct our own plot.
```R
# using DEseq2 built in method
plotMA(deseq2Results)
```

```R
# load libraries
library(ggplot2)
library(scales) # needed for oob parameter

# coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

# plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# let's add some more detail
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)
```
We can see from the above plots that it is in the shape of a trumpet characteristic MA plots. Further we have overlayed density contours in the last plot and as expected these density contours are centered around a 0 ratio. We can further see that as the average counts increase there is more power to call a gene as differentially expressed based on the fold change. You'll also notice that we have quite a few points without an adjusted p value on the left side of the x axis. This is occuring because the results() function automatically performs independent filtering using the mean of normalized counts as a filter statistic. This is done to increase the power to detect an event by not testing those genes which are unlikely to be significant based on their high dispersion.

# Viewing normalized counts for a single geneID
Often it will be usefull to plot the normalized counts for a single gene in order to get an idea of what is occurring at a per sample basis. Fortunately the [plotCounts()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/plotCounts) function from DEseq2 will extract the data we need for plotting.

```R
# extract counts for the gene otop2
otop2Counts <- plotCounts(deseq2Data, gene="ENSG00000183034", intgroup=c("tissueType", "individualID"), returnData=TRUE)

# plot the data using ggplot2
colourPallette <- c("#bbcfc4","#7145cd","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
ggplot(otop2Counts, aes(x=tissueType, y=count, colour=individualID, group=individualID)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1), ) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")
```
From the resulting plot we can see that almost all individuals show down-regulation of this gene in both the primary and met samples compared to the normal. We've also introduced a few new ggplot2 concepts, let's breifly go over them. You will notice that we have specified a [group](http://ggplot2.tidyverse.org/reference/aes_group_order.html) when we initalized our plot. By default ggplot would have assumed the groups were for the discrete variables plotted on the x-axis, when connecting  points with [geom_line()](http://ggplot2.tidyverse.org/reference/geom_path.html) this would have connected all points for each discrete variable instead of connecting by the individual id. Try removing the grouping to get a sense of what happens. We have also altered the legend using [guides()](http://ggplot2.tidyverse.org/reference/guides.html) to specify the legend to act on and [guide_legend()](http://ggplot2.tidyverse.org/reference/guide_legend.html) to specify that the colour legend should have 3 columns for values instead of just 1. Lastly we have added a main title with [ggtitle()](http://ggplot2.tidyverse.org/reference/labs.html). 

### Additional information and references
* [Experimental Data, Kim et al.](https://www.ncbi.nlm.nih.gov/pubmed/25049118)
* [DEseq2](https://www.ncbi.nlm.nih.gov/pubmed/25516281)

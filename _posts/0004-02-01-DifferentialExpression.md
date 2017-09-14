---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Differential expression with DEseq2
categories:
    - Module 4
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-01
---

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Oftentimes, it will be used to define the differences between multiple biological conditions (e.g. drug treated vs. untreated samples). There are many, many tools available to perform this type of analysis. In this course we will rely on a popular Bioconductor package [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). We will then make various visualizations to help interpret our results.

### Dataset
For this analysis we will use the RNAseq data obtained from the (EBI Expression Atlas (GXA))[https://www.ebi.ac.uk/gxa]. Specifically data set [E-GEOD-50760](https://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Downloads) which corresponds to [PMID: 25049118](https://www.ncbi.nlm.nih.gov/pubmed/25049118). This data consists of 54 samples from 18 individuals. Each individual has a primary colorectal cancer sample, a metastatic liver sample, and a normal sample of the surrounding colonic epithilium. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts. We will use the output from HTseq as a starting point. 

It can be downloaded from the GXA.Download these files for use in this exercise:
* Raw counts data from here:[E-GEOD-50760 raw counts](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts)
* Sample information from here: [E-GEOD-50760 sample info](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/ExperimentDesignFile.RnaSeq/experiment-design).

A full description of the experimental design can be found at [array express](http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-50760/) and the [expression atlas](http://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Results).

### How DEseq2 works
[DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a popular differential expression analysis package available through Bioconductor. Its differential expression tests are based on a negative binomial generalized linear model. To get started we will first need to install the package and load the library.
```R
# install the latest version of DEseq2
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
```
### Input data
Input data for [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) consists of non-normalized sequence read counts at either the gene or transcript level. No preliminary normalization of this data is needed. [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) will internally corrects for differences in library size, using the raw counts. The tool [HTseq](http://htseq.readthedocs.io/en/release_0.9.0/) can be used to obtain this information and is what was used for our example data. 

Let's go ahead and load the data and sample information into R. Don't forget to set your working directory to the location of the data files:

```R
# Read in the raw read counts
rawCounts <- read.delim("E-GEOD-50760-raw-counts.tsv")

# Read in the sample mappings
sampleData <- read.delim("E-GEOD-50760-experiment-design.tsv")

# also save a copy for later
sampleData_v2 <- read.delim("E-GEOD-50760-experiment-design.tsv")
```

The next step is to create an object of class DESeqDataSet, which will store the readcounts and intermediate calculations needed for the differential expression analysis. The object will also store the design formula used to estimate dispersion and log2 fold changes used within the model. "Dispersion" is a parameter of the [Generalized Linear Model](https://en.wikipedia.org/wiki/Generalized_linear_model) that relates to to the variance of the distribution. For more details refer to [PMID: 24349066](https://www.ncbi.nlm.nih.gov/pubmed/24349066) and [PMID: 22287627](https://www.ncbi.nlm.nih.gov/pubmed/22287627). 

When specifying the formula it should take the form of a "~" followed by "+" separating factors. When using the default DEseq2 parameters the factor of interest (tissue type in this case) should be specified last and the control within that factor should be first when viewing the [levels()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/levels) for that variable. 

There are 4 methods to create this object depending on the format the input data is in. 

Because we already have our data loaded into R we will use [DESeqDataSetFromMatrix()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class).

```R
# convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

# convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.clinical.information.", "Sample.Characteristic.individual.")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("tissueType", "individualID")
sampleData$individualID <- factor(sampleData$individualID)

# put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$tissueType <- factor(sampleData$tissueType, levels=c("normal-looking surrounding colonic epithelium", "primary colorectal cancer", "metastatic colorectal cancer to the liver"))

# create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ individualID + tissueType)
```

This was quite a bit of code, let's go over whats going on here. The first thing we do is coerce the data frame containing the read counts into a format [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) can accept. Specifically this must be a matrix with row names as genomic features (i.e. genes), and column names as samples. Next [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) requires a data frame specifying the mapping of samples to variables, we load this in and clean it up some keeping only the variables we care about and making sure everything is a factor. For [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to work properly the column names of the count matrix must be in the same order as the row names of the sample mapping data, to ensure this we re-order the column names of the count data and run a check to ensure this has occurred correctly. To take advantage of the default settings of [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) the control of the variable of interest, in our case the tissue type, should be the first element in the levels of that variable. Because we have more than 2 conditions for this variable we will not be taking advantage of the default settings however it's good to get into the practice of doing this so we do it here. We then create a [DEseq2DataSet](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class) object with this information and supply a formula where we use the individual id as a blocking factor and tissue type as the comparison variable.

### Pre-filtering of data
While it is not strictly necessary, it is good to do some preliminary filtering of the data before running the differential expression analysis. This will reduce the size of the [DEseq2DataSet](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class) object and speed up the runtime of the algorithm. Here we are performing relatively minor filtering requiring genes to have more than 5 read of support in any sample.

First see what affect this filter will have
```R
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
```

Now actually apply the filter
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
The next step is to run the function [DEseq()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeq) on our DESeq2 data set object. In this step the algorithm will perform the following:
1. Estimation of size factors
2. Estimation of dispersion
3. Negative Binomial GLM fitting and Wald statistic.

This can take a few minutes to perform, for convenience a .RData object containing the resulting object is available to download [here](http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/deseq2Data_v1.RData). You can load this into your R environment with [load()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/load) either locally after downloading the file or directly through the web.
```R
# run pipeline for differential expression steps
deseq2Data <- DESeq(deseq2Data)

# load the R environment with this object  from the web (optional)
load(url("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/deseq2Data_v1.RData"))

# or directly after downloading the file (optional)
load("deseq2Data_v1.RData")
```

### Extracting results
Finally we can extract the differential expression results with the [results()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results) function. When using this function we need to tell [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) what comparison to make. This is only necessary if the design formula is multi-factorial or, as in our case, the variable in the design formula has more than 2 levels. This is done with the `contrast` parameter which takes a character vector of three elements giving the name of the factor of interest, the numerator (i.e. comparator), and the denominator (i.e. control). Let's get output for normal tissue vs primary tumor  expression results and view a summary of results.
```R
# Extract differential expression results
deseq2Results <- results(deseq2Data, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))

# view summary of results
summary(deseq2Results)
```
### MA-plot
An MA plot plots a log ratio (M) over an average (A) in order to visualize the differences between two groups. In general we would expect the expression of genes to remain consistent between conditions and so the MA plot should be similar to the shape of a trumpet with most points residing on a y intercept of 0. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) has a built in method for constructing an MA plot of our results however as this is a a visualization course let's go ahead and use what we know of [ggplot2](http://ggplot2.tidyverse.org/reference/) to construct our own plots.
```R
# using DEseq2 built in method
plotMA(deseq2Results)
```

```R
# load libraries
# install.packages(c("ggplot2", "scales", "viridis"))
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

# plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# let's add some more detail
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)
```
We can see from the above plots that it is in the shape of a trumpet characteristic MA plots. Further we have overlayed density contours in the last plot and as expected, these density contours are centered around a y-intercept of 0. We can further see that as the average counts increase there is more power to call a gene as differentially expressed based on the fold change. You'll also notice that we have quite a few points without an adjusted p-value on the left side of the x-axis. This is occurring because the [results()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results) function automatically performs independent filtering using the mean of normalized counts. This is done to increase the power to detect an event by not testing those genes which are unlikely to be significant based on their high dispersion.

# Viewing normalized counts for a single geneID
Often it will be useful to plot the normalized counts for a single gene in order to get an idea of what is occurring at a per sample basis. Fortunately the [plotCounts()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/plotCounts) function from DEseq2 will extract the data we need for plotting this.

```R
# extract counts for the gene otop2
otop2Counts <- plotCounts(deseq2Data, gene="ENSG00000183034", intgroup=c("tissueType", "individualID"), returnData=TRUE)

# plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
ggplot(otop2Counts, aes(x=tissueType, y=count, colour=individualID, group=individualID)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")
```
From the resulting plot we can see that almost all individuals show down-regulation of this gene in both the primary tumor and metastasis samples compared to the normal. We've also introduced a few new ggplot2 concepts, so let's briefly go over them. You will notice that we have specified a [group](http://ggplot2.tidyverse.org/reference/aes_group_order.html) when we initialized our plot. By default ggplot would have assumed the groups were for the discrete variables plotted on the x-axis, and when connecting points with [geom_line()](http://ggplot2.tidyverse.org/reference/geom_path.html) this would have connected all points for each discrete variable instead of connecting by the individual id. Try removing the grouping to get a sense of what happens. We have also altered the legend using [guides()](http://ggplot2.tidyverse.org/reference/guides.html) to specify the legend to act on and [guide_legend()](http://ggplot2.tidyverse.org/reference/guide_legend.html) to specify that the colour legend should have 3 columns for values instead of just 1. Lastly we have added a main title with [ggtitle()](http://ggplot2.tidyverse.org/reference/labs.html). You might have noticed that individual 2 looks a little off, specifically based on the plot we created there appears to be two normal samples for this individual. Looking at our sample data we can observe that this is indeed the case, contrary to the stated experimental design. Looking at the sampleData dataframe we created, we can see that our suspicions are correct. Often when visualizing data you will catch potential errors that would have otherwised been missed.

### Visualizing expression data with a heatmap
It is often informative to plot a heatmap of differentially expressed genes and to perform unsupervised clustering based on the underlying data to determine sub categories within the experiment. In this case we can attempt remedy the error we observed in individual 2. We can use [ggplot](http://ggplot2.tidyverse.org/index.html) and [ggdendro](https://cran.r-project.org/web/packages/ggdendro/index.html) for this task, but first we must obtain transformed values from the RNAseq counts. The differential expression analysis started from raw counts and normalized using discrete distributions however when performing clustering we must remove the dependence of the variance on the mean. In other words we must remove the experiment wide trend in the data before clustering. There are two functions within DEseq2 to transform the data in such a manner, the first is to use a regularized logarithm [rlog()](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/rlog) and the second is the variance stablizing transform [vst()](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/vst). There are pros and cons to each method, we will use [vst()](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/vst) here simply because it is much faster. By default both [rlog()](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/rlog) and [vst()](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/vst) are blind to the sample design formula given to [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in [DESeqDataSetFromMatrix()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class). However this is not appropriate if one expects large differences in counts, which can be explained by the differences in the experimental design. In such cases the `blind` parameter should be set to `FALSE`.

```R
# transform count data using the variance stablilizing transform
deseq2VST <- vst(deseq2Data)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)

# Keep only the significantly differentiated genes
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# covert the VST counts to long format for ggplot2
library(reshape2)
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
```

Let's briefly talk about the steps we took to obtain the heatmap we plotted above. First we took our [DESeq2DataSet](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/DESeqDataSet-class) object we obtained from the command [DESeq()](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/DESeq) and transformed the values using the variance stabilizing tranform algorithm from the [vst()](https://www.rdocumentation.org/packages/DESeq/versions/1.24.0/topics/vst) function. We then extracted these transformed values with the [assay()](https://www.rdocumentation.org/packages/SummarizedExperiment/versions/1.2.3/topics/SummarizedExperiment-class) function and converted the resulting object to a data frame with a column for gene id's. We next used the differential expression results we had previously obtained to filter our transformed matrix to only those genes which were significantly differentially expressed (q-value <= .05)and with a log2 fold change greater than 3. Up to this point our transformed values have been in "wide" format, however [ggplot2](http://ggplot2.tidyverse.org/reference/) required long format, we achieve this with [melt()](https://www.rdocumentation.org/packages/reshape2/versions/1.4.2/topics/melt) function from the [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html) package. Finally we use [ggplot2](http://ggplot2.tidyverse.org/reference/) and the [geom_raster()](http://ggplot2.tidyverse.org/reference/geom_tile.html) function to create a heatmap using the color scheme available from the [viridis](https://cran.r-project.org/web/packages/viridis/index.html) package.

Now that we have a heatmap let's start clustering using the functions available with base R. The first step is to convert our transformed values back to wide format using [dcast()](https://www.rdocumentation.org/packages/reshape2/versions/1.4.2/topics/cast) with row names and columns names as genes and samples respectively. Next we can compute a distance matrix using [dist()](https://www.rdocumentation.org/packages/stats/versions/3.4.1/topics/dist) which will compute the distance based on the rows of our data frame. This will give us the distance matrix based on our genes, but we also want a distance matrix based on the samples, for this we simply have to transpose our matrix before calling [dist()](https://www.rdocumentation.org/packages/stats/versions/3.4.1/topics/dist) with the function [t()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/t). From there we can perform hierarchical clustering with the [hclust()](https://www.rdocumentation.org/packages/stats/versions/3.4.1/topics/hclust) function. We then install the [ggdendro](https://cran.r-project.org/web/packages/ggdendro/index.html) package to construct a dendrogram as described in their [vignette](https://cran.r-project.org/web/packages/ggdendro/vignettes/ggdendro.html). Finally we re-create the heatmap to match the dendrogram using the [factor()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/factor) function to re-order how samples are plotted. Because ggplot plots use grid graphics underneath we can use the [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) package to combine both plots into one with the function [grid.arrange()](https://www.rdocumentation.org/packages/gridExtra/versions/2.2.1/topics/arrangeGrob).

```R
# convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# compute a distance calculation on the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# construct a dendogram for samples
install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# construct the heatmp
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())

# combine the dendrogram and the heatmap
install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))
```

Our graph is looking pretty good, but you'll notice that the two plots don't seem to line up. This is because the plot widths from our two plots don't quite match up. This can occur for a variety of reasons however in this case it is because we have a legend in one plot but not in the other. Fortuanately this sort of problem is generally easy to fix.

```R
# load in libraries necessary for modifying plots
#install.packages("gtable")
library(gtable)
library(grid)

# modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

# add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# draw the plot
grid.draw(finalGrob)
```

You will notice that we are loading in a few additional packages to help solve this problem. All plots within [ggplot](http://ggplot2.tidyverse.org/reference) are at the lowest level graphical objects (grobs) from the [grid]() package. The [gtable](https://cran.r-project.org/web/packages/gtable/index.html) allows us to view and manipulate these grobs as tables making them easier to work with. The first step after loading these packages is to alter the x-axis scales in our plots. By default [ggplot](http://ggplot2.tidyverse.org/reference) adds a padding within all plots so data is not plotted on the edge of the plot. We can control this padding with the `expand` parameter within ggplot's [scale](http://ggplot2.tidyverse.org/reference/scale_continuous.html) layers. To get things just right you will need to alter these parameters so data within the plots line up, aproximate alterations are provided above so let's move on the the next step. We need to make sure the main panels in both plots line up exactly, in order to do this we must first convert these plots to table grobs to make them easier to manipulate, this is achieved with the [ggplotGrob()](https://www.rdocumentation.org/packages/ggplot2/versions/2.2.1/topics/ggplotGrob) function. Now that these are grobs we compare the widths of each grob, you might notice that two widths seem to be missing from the the dendrogram plot. These two widths correspond to the legend in the heatmap, we use the [gtable_add_cols()](https://www.rdocumentation.org/packages/gtable/versions/0.2.0/topics/gtable_add_cols) function to insert in the widths for these two elements from the heatmap grob into the dendrogram grob. The next step is to then find the maximum widths for each grob object using [unit.pmax()](https://www.rdocumentation.org/packages/grid/versions/3.4.1/topics/unit.pmin) and to overwrite each grobs width with these maximum widths. Finally we arrange the grobs on the page we are plotting to using [arrangeGrob()](https://www.rdocumentation.org/packages/gridExtra/versions/2.2.1/topics/arrangeGrob) and draw the final plot with [grid.draw()](https://www.rdocumentation.org/packages/grid/versions/3.4.1/topics/grid.draw).

Now that we've completed that let's add a plot between the dendrogram and the heatmap showing the tissue type to try and discern which normal sample in individual 2 should be the metastasis sample as per the experimental design described.

```R
# re-order the sample data to match the clustering we did
sampleData_v2$Run <- factor(sampleData_v2$Run, levels=clusterSample$labels[clusterSample$order])

# construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52")
sampleClinical <- ggplot(sampleData_v2, aes(x=Run, y=1, fill=Factor.Value.clinical.information.)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tissue", values=colours) + theme_void()

# convert the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

# make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)
```

{% include question.html question="Based on the plot you produced above, which normal sample for individual 2 is likely the metastasis sample" answer='Based on the clustering, sample SRR975587'%}

### Exercises

In the tutorial above we have completed the steps necessary to produce a dendrogram for samples, we could also do this for genes. Try to reproduce the plot below, if you get stuck you'll find an rscript that below which will produce the correct answer! Your steps should look something like this:

{% include figure.html image="/assets/Deseq2/heatmap_deseq2_final.png" position="center" %}

1. Create a dendrogram for genes
    1. use [coord_flip()](http://ggplot2.tidyverse.org/reference/coord_flip.html) and [scale_y_reverse()](http://ggplot2.tidyverse.org/reference/scale_continuous.html)
2. re-arrange the gene cells on the y-axis to match the gene dendrogram
3. convert the dendrogram and heatmaps to grobs
4. align the gene dendrogram heights to match the heatmap (note: you don't need to worry about the legend in this case)
    1. you'll also need to repeat the tutorial steps for aligning the sample dendrogram since we have a new heatmap, consider this step 4B.
5. create a blank panel you can use this: `blankPanel<-grid.rect(gp=gpar(col="white"))`
6. use the gridExtra package to position the panels and plot the result.

{% include question.html question="Reproduce the plot above, follow the steps outlined." answer='This Rscript <a href="http://genomedata.org/gen-viz-workshop/intro_to_deseq2/exercise_1/exercise1_deseq2_final_heatmap.R">file</a> contains the correct answer.'%}

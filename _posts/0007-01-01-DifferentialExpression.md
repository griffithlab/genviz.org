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
While not strictly necessary it is good to do some preliminary filtering of the data before running the differential expression analysis. This will reduce the size of the DESeq2 object and speed up the the speed of the algorithm. Here we are performing relatively minor filtering requiring genes to have more than 1 read of support.
```R
# perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 1, ]
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

This step can take a few minutes to perform, for convenience a .RData object containing an R environment up to this step is available to download here.
```R
deseq2Data <- DESeq(deseq2Data)
```

### Extracting results
Finally we can extract the differential expression results with the [results()](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results) function. When using this function we need to tell DESeq2 what comparison to make. This is only neccessary if the design formula is multi-factorial or as in our case the variable in the design formula has more than 2 levels. This is done with the contrast parameter which takes a character vector of three elements giving the name of the factor of interest, the numerator, and the denominator (i.e. control). Let's get output for normal vs primary expression results.
```R
# Extract differential expression results
deseq2Results <- results(deseq2Data, contrast=c("tissueType", primary colorectal cancer", "normal-looking surrounding colonic epithelium"))
```

### Additional information and references
* [Experimental Data, Kim et al.](https://www.ncbi.nlm.nih.gov/pubmed/25049118)
* [DEseq2](https://www.ncbi.nlm.nih.gov/pubmed/25516281)

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
```
### Input data
Input data for [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) consists of un-normalized sequence read counts at either the gene or transcript level. No normalization of this data is needed because [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) internally corrects for factors, specifically library size using these raw counts. The tool [HTseq](http://htseq.readthedocs.io/en/release_0.9.0/) can be used to obtain this information and is what was used for our [example data](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts), let's go ahead and load this data and the [sample information](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/ExperimentDesignFile.RnaSeq/experiment-design) into R.

```R
# Read in the raw read counts
rawCounts <- read.delim("E-GEOD-50760-raw-counts.tsv")

# Read in the sample mappings
sampleData <- read.delim("E-GEOD-50760-experiment-design.tsv")
```

### Additional information and references
* [Experimental Data, Kim et al.](https://www.ncbi.nlm.nih.gov/pubmed/25049118)
* [DEseq2](https://www.ncbi.nlm.nih.gov/pubmed/25516281)

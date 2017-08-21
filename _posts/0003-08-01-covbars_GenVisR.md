---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to Sequencing Coverage Plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-08-01
---

Commonly when a sample has undergone sequencing you will want to know the sequencing depth achieved in order to get an idea of the data quality. The [covBars()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/covBars) function from the GenVisR package is designed to help in visualizing this sort of data by constructing a color ramp of cumulative coverage. In this section we will be reconstructing the capture coverage plots for 4 samples from the paper ["Comprehensive genomic analysis reveals FLT3 activation and a therapeutic strategy for a patient with relapsed adult B-lymphoblastic leukemia"](https://www.ncbi.nlm.nih.gov/pubmed/27181063) shown in Supplemental Figure S3.

{% include figure.html image="/assets/GenVisR/Coverage_Summary.png" width="650" link="http://www.sciencedirect.com/science/article/pii/S0301472X16301151?via%3Dihub#appsec2" title="Sequence Coverage" author="Griffith et al." license="CC BY-NC-ND" license_link="https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode" %}

### Data Preparation
The [covBars()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/covBars) function takes as input a matrix with rows representing the sequencing depth, columns representing samples, and matrix values representing the number of reads meeting that criteria for a given cell. The best way to begin constructing this data is with the command line tool [samtools depth](http://www.htslib.org/doc/samtools.html) for anything other than whole genome sequencing data in which case [picards CollectWgsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html) program might be a better choice. The output from [samtools depth](http://www.htslib.org/doc/samtools.html) for 4 capture samples from the adult B-lymphoblastic leukemia manuscript linked above is available to download [here](http://genomedata.org/gen-viz-workshop/GenVisR/ALL1_CaptureDepth.tsv). The command used to create this file was `samtools depth -q 20 -Q 20 -b roi.bed -d 12000 sample1.bam sample2.bam sample3.bam sample4.bam`, let's briefly go over the parameters used to create this file and load it into R. Within the [samtools depth](http://www.htslib.org/doc/samtools.html) command we used the parameters `-q 20` and `-Q 20` to set a minimum base and mapping quality respectively. Essentially this means that a read will not be counted if these criteria are not met. We used `-b roi.bed` to tell [samtools depth](http://www.htslib.org/doc/samtools.html) to only create pileups for those  regions encompassed by each line of the bed file. We used `-d 12000` to ensure pileups are not stopped until 12,000 reads are reached for that coordinate. Finally we specified the four bam files for which we want to create pileups for, each bam file has been indexed with [samtools index](http://www.htslib.org/doc/samtools.html).

```R
# load the ouput of samtools depth into R
seqData <- read.delim("~/Desktop/ALL1_CaptureDepth.tsv", header=F)
colnames(seqData) <- c("chr", "position", "samp1", "samp2", "samp3", "samp4")
seqData <- seqData[,c(3:6)]
```

### Manipulating data
After reading in the data you'll notice that the format is not quite what [covBars()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/covBars) will accept; Instead of having read pileup counts for each coverage value (1X 2X 3X etc.) we have read pileups for each coordinate position in the bed file. Let's go step by step how to manipulate the data into the proper format for [covBars()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/covBars). The first thing we need to do loop over all columns of `seqData` and obtain a tally of read pileups for each coverage value, this can be done with a call to [apply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/apply) using the [plyr](https://cran.r-project.org/web/packages/plyr/index.html) function [count()](https://www.rdocumentation.org/packages/plyr/versions/1.8.4/topics/count). Strictly for readability we then create our own function called `renameCol` with the purpose of renaming the columns names to a more human readable format and use [lapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply) to apply the function to each data frame in our list. You might notice that each data frame in `seqCovList` has a differing number of rows. This has occurred because not every sample will have pileups for the same coverage values, this is especially true when getting into higher coverage depths owing to outliers. We will fix this by creating a framework data frame containing all possible coverage values, from the min to the max observed from each data frame in the list `seqCovList`. For this we use [lapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply) to apply [max()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extremes) to each coverage column within the data frames in `seqCovList` using an anonymous function within the [lapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply). This will return a list of the maximum coverage value for each data frame so we use [unlist()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/unlist) to coerce the list to a vector and use [max()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extremes) again to find the overall maximum coverage observed. A similar procedure is perormed to obtain the minimum coverage and we create are framework data frame using these maximum and minimum values with [data.frame()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/data.frame). Now that we have our framework data frame we can use [merge()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/merge)  and [lapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply) to perform a natural join between the the framework data frame `covFramework` and each data frame containg the acutal read pileups for each coverage in `seqCovList`.

For this we use the functions [min](), [max](), [lapply](), and [unlist]().
```R
# Count the occurrences of each coverage value
# install.packages("plyr")
library(plyr)
seqCovList <- apply(seqData, 2, count)

# rename the columns for each dataframe
renameCol <- function(x){
    colnames(x) <- c("coverage", "freq")
    return(x)
}
seqCovList <- lapply(seqCovList, renameCol)

# create framework data frame with entries for the min to the max coverage
maximum <- max(unlist(lapply(seqCovList, function(x) max(x$coverage))))
minimum <- min(unlist(lapply(seqCovList, function(x) min(x$coverage))))
covFramework <- data.frame("coverage"=minimum:maximum)

# Merge the framework data frame with the coverage
# seqCovList <- lapply(seqCovList, function(x, y) merge(x, y, by="coverage", all=TRUE), covFramework)
seqCovList <- lapply(seqCovList, merge, covFramework, by="coverage", all=TRUE)

# merge all data frames together
seqCovDataframe <- Reduce(function(...) merge(..., by="coverage", all=T), seqCovList)

# set all NA values to 0
seqCovDataframe[is.na(seqCovDataframe)] <- 0

# set the rownames, remove the extra column, and convert to a matrix
rownames(seqCovDataframe) <- seqCovDataframe$coverage
seqCovDataframe$coverage <- NULL
seqCovMatrix <- as.matrix(seqCovDataframe)

# rename columns
colnames(seqCovMatrix) <- c("sample1", "sample2", "sample3", "sample4")
```

```R
# run covbars
covBars(seqCovMatrix)
```

---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to GenVisR
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-01-01
---

The advent of next generation sequencing (NGS) has allowed for the production of massive amounts of genomic data that is available for analysis. While there exist bioinformatics tools to identify mutational events found within individual samples, these tools are insufficient to fully understand the biological impact these variants have in a disease like cancer. To help alleviate this bottleneck, GenVisR is a bioconductor package in R that allows for the generation of highly-customized publication quality visualizations that aid in the analysis of cohort-level genomic data. In this module, we will become familiar with the GenVisR package and work with several functions to visualize/analyze genomic variations in the form of single nucleotide variants (SNVs), insertions, deletions, loss of heterozygosity, and copy number variants. We will also look into a function that can visualize the coverage sequencing depth acquired across all samples in a cohort. For each function, we will work with the fundamental parameters needed to generate plots that convey enough information for analysis. A full list of parameters with their description, and other functions as part of the GenVisR package can be found on the GenVisR vignette from the Bioconductor website.
 [https://bioconductor.org/packages/release/bioc/html/GenVisR.html](https://bioconductor.org/packages/release/bioc/html/GenVisR.html).

```R
# install the GenVisR package from Bioconductor. Try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GenVisR")

# load the package into R
library(GenVisR)
```

## More GenVisR Functions
For a full list of GenVisR functions with descriptions and examples, see the GenVisR reference manual and vignette on the Bioconductor page at [https://bioconductor.org/packages/release/bioc/html/GenVisR.html](https://bioconductor.org/packages/release/bioc/html/GenVisR.html).

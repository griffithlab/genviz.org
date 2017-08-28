---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to GenVisR
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-01-01
---

The advent of next generation sequencing (NGS) has allowed for the production of massive amounts of genomic data that is available for analysis. There are many tools available for the analysis and visualization of this type of data. In this module we will focus on the R package [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html), The reasons for doing so can be summarized as follows:
1. The package is built upon [ggplot2](http://ggplot2.tidyverse.org/reference/) and will allow us to leverage information we've learned in previous modules.
2. The package is intended to be flexible supporting multiple file types, species, etc.
3. The package is relatively popular; in the top 20% of biooconductor downloads.
4. The package is regularly updated with improvements, bug fixes, etc.
5. The package is maintained by the griffithlab

We reccomend installing [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) from bioconductor in order to ensure the most stable version.

```R
# install the GenVisR package from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("GenVisR")

# load the package into R
library(GenVisR)
```

### More GenVisR Functions
For a full list of GenVisR functions with descriptions and examples, see the GenVisR reference manual and vignette on the Bioconductor page at [https://bioconductor.org/packages/release/bioc/html/GenVisR.html](https://bioconductor.org/packages/release/bioc/html/GenVisR.html).

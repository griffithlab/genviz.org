---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to GenVisR
categories:
    - Module-03-GenVisR
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-01-01
---

The advent of next generation sequencing (NGS) has allowed for the production of massive amounts of genomic data. There are many methods and tools available for the analysis and visualization of these data. In this module we will focus on the R package [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html). GenVisR offer several advantages for visualization and interpretation of genomic data:
1. GenVisR is built upon [ggplot2](http://ggplot2.tidyverse.org/reference/) and thus allows the user to leverage the many existing graphical functions of that package (as well as the information we've learned in previous modules).
2. GenVisR is intended to be flexible, supporting multiple common genomic file formats, species, etc.
3. GenVisR attempts to make popular, but very complex genomic visualizations much simpler to produce. It essentially offers convenience functions that allow sophisticated plots to be made from common genomic data types with just a handful of lines of code.
4. GenVisR is relatively popular (in the top 20% of bioconductor downloads) and therefore benefits from a large community of users, many published examples, and a number of online tutorials.
5. GenVisR is maintained by the [griffithlab](https://github.com/griffithlab/GenVisR) and is regularly updated with improvements, bug fixes, etc.

We recommend installing [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) from bioconductor in order to ensure the most stable version.

```R
# Install GenVisR
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenVisR", version = "3.8")

#NOTE: If you have an older version of R/Bioconductor you may need to use the previous Bioconductor install method
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenVisR")

# load the package into R
library(GenVisR)
```

### More GenVisR Functions
For a full list of GenVisR functions with descriptions and examples, see the GenVisR reference manual and vignette on the Bioconductor page at [https://bioconductor.org/packages/release/bioc/html/GenVisR.html](https://bioconductor.org/packages/release/bioc/html/GenVisR.html).

### Getting Help

* Visit [GenVisR page at BioConductor](https://bioconductor.org/packages/release/bioc/html/GenVisR.html)
* Read [GenVisR Reference Manual](https://bioconductor.org/packages/release/bioc/manuals/GenVisR/man/GenVisR.pdf)
* In R type `browseVignettes` to view available vignettes: [Introduction](https://bioconductor.org/packages/release/bioc/vignettes/GenVisR/inst/doc/Intro.html), [Visualizing Small Variants](https://bioconductor.org/packages/release/bioc/vignettes/GenVisR/inst/doc/Upcoming_Features.html), and [Waterfall](https://bioconductor.org/packages/release/bioc/vignettes/GenVisR/inst/doc/waterfall_introduction.html)
* BioStar Tutorials: [waterfall](https://www.biostars.org/p/181159/) and [lolliplot](https://www.biostars.org/p/196576/) 

[Module 3 Lecture](https://github.com/griffithlab/gen-viz-lectures/raw/master/GenViz_Module3_Lecture.pdf)

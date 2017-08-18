---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to Copy Number Spectrum Plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-05-01
---

## cnSpec
Large-scale chromosomal rearrangements can occur in the form of copy number variants (CNV) with loss of function (DEL) or gain of function (AMP) consequence. The cnSpec function plots regions of CNV loss and amplification across the genome of all samples in a cohort. This allows for the quick identification of regions that are recurrently amplified or deleted. A dataframe consisting of the column names and values for chromosome, start, end, segmean, and sample serves as the necessary input for this function. The segmean variable contains values that indicate regions of CNV (can be obtained through CBS algorithm and DNAcopy). To specify the boundaries for each chromosome, the UCSC genome 'hg19' is used by default. When working with non-human genomes, or with UCSC genome hg38, other genomes can be specified. These include hg38, mm10, mm9, and rn5. If the UCSC genome assembly is still not one of these, cnSpec will query the UCSC sql database to obtain the necessary chromosomal boundary information. Alternatively, the chromosomal boundaries can be specified by the parameter, y. Input for this argument is a dataframe containing the columns: chromosome, start, and end.

<!--Example taken from GenVisR vignette-->
```R
#hg19chr dataframe defines the choromosomal boundaries
head(hg19chr)

# Generate cnSpec plot
cnSpec(lucCNseg, y=hg19chr)
```

<!--Talk about difference between absolute and relative-->
As seen in the above plot, CN values are color-coded based on a gradient where segment values below 2 are blue (indicating a deletion) and segment values above 2 are red (indicating an amplification). This represents an absolute scale where copy neutral events are specified by a segment value of 2 (there are 2 alleles). Sometimes, CN data is provided on a relative scale, where copy neutral events correspond to a segment value of 0, where negative values represent a segment of copy number loss and positive values represent a segment of copy number gain. The absolute and relative scale can be specified with the <font color="green">CNscale</font> parameter.

<!--Example taken from GenVisR vignette-->
```R
# Generate a cnSpec plot on a relative scale
cnSpec(lucCNseg, CNscale = "relative")
```

---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to gene coverage plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-04-01
---

Often is is usefull to view coverage of a specific region of the genome in the context of specific samples. During initial stages of analysis this can be done with a genome browser such as [IGV](http://software.broadinstitute.org/software/igv/) however when preparing a publication more fine grain control is usefull. For example you may wish to change the coverage scale, reduce the size of introns, or visualize many samples at once. GenVisR provides a function for this aptly named [genCov](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/genCov).

### introduction to datasets

In this section we will be using coverage data derived from two mouse samples from the study ["Truncating Prolactin Receptor Mutations Promote Tumor Growth in Murine Estrogen Receptor-Alpha Mammary Carcinomas"](https://www.ncbi.nlm.nih.gov/pubmed/27681435). We will be showing that the knockout of the *STAT1*
 gene described in the manuscript was successful. In order to obtain the preliminary data we used the command `bedtools multicov -bams M_CA-TAC245-TAC245_MEC.prod-refalign.bam M_CA-TAC265-TAC265_MEC.prod-refalign -bed stat1.bed` to obtain coverage values for the wildtype TAC245 sample and the tumor free knockout TAC265 sample.

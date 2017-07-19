---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Differential Expression with DEseq2
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0007-01-01
---

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Oftentimes, it will be used to define the differences between multiple biological conditions (e.g. post-treatment vs. untreated samples). There are many tools available to perform this type of analysis, in this course we will rely on the bioconductor package (DEseq2)[https://bioconductor.org/packages/release/bioc/html/DESeq2.html]. We will then make various visualizations to help interpret our results.

### dataset
For this analysis we will use the RNAseq data from (E-GEOD-50760)[https://www.ncbi.nlm.nih.gov/pubmed/25049118]. This data consists of 54 samples from 18 individuals, each individual has a primary colorectal cancer sample, a metastatic liver sample, and a normal sample of the surrounding colonic epithilium. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts. We will use the output from HTseq as a starting point which can be downloaded (here)[https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts], sample information can be downloaded (here)[https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/ExperimentDesignFile.RnaSeq/experiment-design]. A full description of the experimental design can be found at (array express)[http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-50760/] and the (expression atlas)[http://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Results].

### dEseq2

Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.


### Additional information and references
* [Experimental Data, Kim et al.](https://www.ncbi.nlm.nih.gov/pubmed/25049118)

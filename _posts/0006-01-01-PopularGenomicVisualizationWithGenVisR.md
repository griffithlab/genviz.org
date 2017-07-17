---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to GenVisR
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-01-01
---

<!-- Introduce GenVisR and the problem it solves -->
Over the years, the efficiency of __. This has generated a massive amount of data (on the level of __ ) that is availabe for analysis. To date, there are insufficient bioinformatics tools that can __. To alleviate this bottleneck, GenVisR is a bioconductor package in R that aids in the analysis of cohort-level genomic data. In this module, we will become familiar with the GenVisR package and work with several functions to visualize genomic data at the level of single nucleotide variants (SNVs), insertions, deletions, loss of heterozygosity, copy number variants, and coverage. 

The GenVisR package from the Bioconductor website and following the installation instructions: [https://bioconductor.org/packages/release/bioc/html/GenVisR.html] (https://bioconductor.org/packages/release/bioc/html/GenVisR.html).

<!-- Install GenVisR in Rstudio 
```R
# install the GenVisR package from Bioconductor. Try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GenVisR")
```
<!-- Talk about waterfall function -->
### waterfall
In all genomic analysis projects, it is important for a researcher to know which genes are frequently mutated in the cohort and have an idea on the biological consequence of each mutation. The waterfall function allows for the vislualization of this mutational landscape from all samples in a cohort. As seen on the right, the main panel portrays the mutation occurence and type, which is indicated by the color of the filled rectangle. Subplots illustrate the percentage of samples with each mutatiion, the mutaitonal burden in each sample, and clinical variables that describe each sample. Often, conflicts arise where multiple mutations in the same gene/sample cell are reported by variant callers and subsequently provided into the input file. These conflicts are resolved by identiyfing the most deleterious mutation as defined by the order of the "mutation type" legend.



Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### tvti
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### cnSpec
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### cnView
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### lohSpec
lohSpec uses a sliding window approach to calcualte the average variant allele frequency (VAF) difference between matched tumor normal samples. This calculation is repeated across the entire genome of all samples in a cohort. 

### genCov
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

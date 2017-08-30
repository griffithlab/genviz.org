---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to Loss of Heterozygosity Plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-07-01
---

We've gone through visualizations of point mutations and copy number changes using [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html), another type of genomic alteration that is often usefull to visualize is ["Loss of Heterozygosity"](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) events. In a diploid genome there are pairs of chromosomes each containing a single copy of the genome. These pairs each come from haploid gametes and are slightly different from each other leading to heterozygosity throughout much of the genome. Situations can arise however where this inherit heterozygosity of the genome is loss, commonly this is through deletions of a parental copy within a chromosome region also known as hemizygousity. Viewing deletions however will not give a complete picture of [LOH](https://en.wikipedia.org/wiki/Loss_of_heterozygosity) as events can arise that will lead to copy-neutral LOH, for example if a parental copy was deleted but the then remaining copy underwent an amplification. In this section we will use the [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) function [lohSpec](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/lohSpec) created specifically for viewing loh events.

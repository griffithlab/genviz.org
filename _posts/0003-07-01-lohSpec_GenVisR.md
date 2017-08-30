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
<<<<<<< HEAD

lohSpec uses a sliding window approach to calculate the average variant allele frequency (VAF) difference between matched tumor normal samples. Essentially, a window of a specified size (defined by <font color="green">window_size</font>) will obtain the absolute mean difference between the tumor and normal VAF (set at 0.5) for each genomic position within the window's position. This window will then move forward by a length specified by the parameter <font color="green">step</font>, and again cauclate the absolute mean tumor-normal VAF difference. This value across all overlapping windows is averaged. Each value is then stored for the overlapping region and this process is repeated for each chromosome in all samples in the cohort. Plotting the X and Y chromosome is optional but can be included when a character vector containing the gender of each sample ("M" for male and "F" for female) is supplied to the <font color="green">gender</font> parameter.

We will be using lohSpec to visualize the landscape of LOH in breast cancer?  The data to generate this figure can be found [here](http://genomedata.org/gen-viz-workshop/GenVisR/HCC1395.varscan.tsv). This dataset is the output from a bioinformatics tool called [varscan()](http://varscan.sourceforge.net/) that can broadly detect mutational events in the form of germline, somatic, and copy number alterations.

This data however contains a lot of unnecessary information. Knowing what you now know about LOH, can you identify the unnecessary pieces of information in the file? 
{% answer="The file contains somatic calls, which do not represent LOH events." %}

To separate out this dataset and obtain the rows that are necessary for our analysis, please open up your terminal (if you are on a mac) or __ (if you are on windows). Enter the following command:
```R
grep "LOH\|Germline" ~/YOURDIRECTORY/HCC1395.varscan.snvs.hq > HCC1395.varscan.germline.loh.tsv
```
A few more steps need to be done to this dataset so that it is usable by lohSpec. As it currently stands, lohSpec requires the user to convert all percentages to decimals (i.e. 33.33% must be 0.33). Additionally, lohSpec requires the dataset to contain specific column names. Fortunately, this can all be done in R. 

```R
## First read in the data
loh_data <- read.delim("HCC1395.varscan.germline.loh.tsv")

## Assign column names that are required 
colnames(loh_data) <- c("chromosome", "position", "n_vaf", "t_vaf", "sample")

## Add the sample name and change the vaf percentages to decimal values
loh_data$sample <- "HCC1395"
loh_data$n_vaf <- as.numeric(gsub("%", "", loh_data$n_vaf))/100
loh_data$t_vaf <- as.numeric(gsub("%", "", loh_data$t_vaf))/100

## Fix the chromosome values
temp$chromosome <- gsub(" ", "", as.character(temp$chromosome))

## Try to run lohSpec
lohSpec(x=loh_data)

The output should report an error message along the lines of: "Detected values with a variant allele fraction either above .6 or below .4 in the normal. Please ensure variants supplied are heterozygous in the normal!" 
```
Given what you know about loh, does this make sense? Variants with a VAF value below 0.4 and above 0.6 in the normal repreesnt germline LOH events. In other words, the normal HCC1395 sample contains LOH at these genomic positions What we are trying to do is plot the somatic LOH landscape observed in the HCC1395 sample, so it is important to remove these variants from the dataset.

```R
## Obtain variants with a VAF greater than 0.4 and below 0.6. 
loh_data <_ subset(loh_data, loh_data$n_vaf > 0.4 & loh_data$n_vaf < 0.6)

## Get autosomes and sex chromosomes (remove GL contigs and mitochondrial chromosomes)
chr <- c(seq(1:22, "X", "Y"))
loh_data <- loh_data[which(loh_data$chromosome %in% chr),]

## Now run lohSpec 
lohSpec(x=loh_data)
``` 
{% include figurehtml image="/assets/GenVisR/lohSpec.vq.png" width="750" %}

<!--Talk about caveats with LOH-->
When interpreting this visualization and other LOH data, it is important to understand that not all LOH corresponds to an actual deletion event. For instance, an amplification of a specific allele/genomic region can appear as LOH since one region will have significantly more copies than the corresponding region on the other chromosome. Additionally, copy neutral loss of heterozygosity can occur where an individual receives a chromosomal pair from one parent and nothing from the other parent. While this situation may appear as LOH since there is only a single chromosomal type present, the individual still has 2 chromosomes (hence neutral copy number). Therefore, it is important to look at coverage to see LOH truly results in a deleted region.
=======
>>>>>>> a2852c474ad3826698beb0bcfbce4340bea14c1f
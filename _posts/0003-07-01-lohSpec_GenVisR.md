---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to Loss of Heterozygosity Plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-07-01
---

## lohSpec (IN THE PROCESS OF MAKING THIS OBJECT ORIENTED)
lohSpec uses a sliding window approach to calcualte the average variant allele frequency (VAF) difference between matched tumor normal samples. Essentially, a window of a specified size (defined by <font color="green">window_size</font>) will obtain the absolute mean difference between the tumor and normal VAF (set at 0.5) for each genomic position within the window's position. This window will then move forward by, specified by the parameter <font color="green">step</font>, and again cauclate the absolute mean tumor-normal VAF difference. This value across all overlapping windows is averaged. Each value is then stored for the overlapping region and this process is repeated for each chromosome in all samples in the cohort. Plotting the X and Y chromosome is optional but can be included when a character vector containing the gender of each sample ("M" for male and "F" for female) is supplied to the <font color="green">gender</font> parameter.

<!--Example taken from GenVisR vignette-->
```R
## Generate LOHspec with default parameters
lohSpec(x=HCC1395_Germline)
````

<!--Talk about caveats with LOH-->
When interpreting this visualization and other LOH data, it is important to understand that not all LOH corresponds to an actual deletion event. For instance, an amplification of a specific allele/genomic region can appear as LOH since one region will have significantly more copies than the corresponding region on the other chromosome. Additionally, copy neutral loss of heterozygosity can occur where an individual receives a chromosomal pair from one parent and nothing from the other parent. While this situation may appear as LOH since there is only a single chromosomal type present, the individual still has 2 chromosomes (hence neutral copy number). Therefore, it is important to look at coverage to see LOH truly results in a deleted region.

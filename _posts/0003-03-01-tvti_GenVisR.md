---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to transition/transverion plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-01
---

## tvti
To better characterize single nucleotide variants (SNVs), it is important to know the frequency of transitions and transversions across the cohort. The TvTi function in GenVisR uses a stacked bar plot to illustrate the proportion/frequency of transition and transversions (specified by the <font color = "green">fileType</font> argument). Similar to the waterfall plot, the input to the tvti function is a maf file or a custom input dataframe with the column names: "sample", "reference", and "variant". Multinucleotide substitutions will be excluded from the plot. To change the color of the tvti bars, change the <font color = "green">pallete</font> parameter and supply it with a character vector with 6 different colors to represent each of the transition and tranversion types.

<!--Example taken from GenVisR vignette-->
```R
## Generate TvTi plot and visualize Proportion/Frequency of each transition and transversion
tvti_type <- "Proportion"
TvTi(brcaMAF, lab_txtAngle=75, fileType="MAF, type = tvti_type)
```

What makes TvTi even more powerful as a genomic visualization tool is the option for the user to plot the expected frequency/proportion of transitions and transversions, and compare the cohort data to these expected values. The expected rate of transitions and transversions is supplied to the parameter <font color = "green">y</font> as a vector with names that match each transition and transversion type. These names must be:  “A->C or T->G (TV)”, “A->G or T->C (TI)”, “A->T or T->A (TV)”, “G->A or C->T (TI)”, “G->C or C->G (TV)”, and “G->T or C->A (TV)”. If this vector is supplied to <font color = "green">y</font>, an additional subplot will be generated to illustrate the expected transition and transversion rates.

<!--Example taken from GenVisR vignette-->
```R
## Fill out vector with expected rates for each transition and transversion type
expec <- c(`A->C or T->G (TV)` = 0.055, `A->G or T->C (TI)` = 0.217, `A->T or T->A (TV)` = 0.076,
    `G->A or C->T (TI)` = 0.4945, `G->C or C->G (TV)` = 0.0645, `G->T or C->A (TV)` = 0.093)

## Illustrate expected transition and transversion rates
TvTi(brcaMAF, y=expec, lab_txtAngle = 45, fileType = "MAF")
```

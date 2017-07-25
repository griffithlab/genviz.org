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
The advent of next generation sequencing (NGS) has allowed for the production of massive amounts of genomic data that is available for analysis. While there exist bioinformatics tools to identify mutational events found within individual samples, these tools are insufficient to fully understand the biological impact these variants have in a disease like cancer. To help alleviate this bottleneck, GenVisR is a bioconductor package in R that allows for the generation of highly-customized publication quality visualizations that aid in the analysis of cohort-level genomic data. In this module, we will become familiar with the GenVisR package and work with several functions to visualize/analyze genomic variations in the form of single nucleotide variants (SNVs), insertions, deletions, loss of heterozygosity, and copy number variants. We will also look into a function that can visualize the coverage sequencing depth acquired across all samples in a cohort. For each function, we will work with the fundamental parameters needed to generate plots that convey enough information for analysis. A full list of parameters with their description, and other functions as part of the GenVisR package can be found on the GenVisR vignette from the Bioconductor website. [https://bioconductor.org/packages/release/bioc/html/GenVisR.html](https://bioconductor.org/packages/release/bioc/html/GenVisR.html).
<!-- Install GenVisR in Rstudio -->
```R
# install the GenVisR package from Bioconductor. Try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GenVisR")

# load the package into R
library(GenVisR)
```
<!-- Talk about waterfall function -->
## Waterfall Function
When conducting genomic analysis, it is important for a researcher to know which genes are frequently mutated in the cohort and have an idea on the biological consequence of each mutation. The waterfall function allows for the vislualization of this mutational landscape from all samples in a cohort. As seen on the right, the main panel portrays the mutation occurence and type, which is indicated by the color of the filled rectangle. Subplots illustrate the percentage of samples with each mutatiion, the mutaitonal burden in each sample, and clinical variables that describe each sample. Often, conflicts arise where multiple mutations in the same gene/sample cell are reported by variant callers and subsequently provided into the input file. These conflicts are resolved by identiyfing the most deleterious mutation as defined by the order of the "mutation type" legend. 

#### Waterfall - Input Data
The input data for the waterfall function can come from a .maf file, which was discussed in Module 1: Introduction to Genomic Visualization and Interpretation ([http://genviz.org/module%201/0001/01/01/Introduction_to_Genomic_Visualization_and_interpretation/](http://genviz.org/module%201/0001/01/01/Introduction_to_Genomic_Visualization_and_interpretation/)). The MAF filetype must be specified in the waterfall function with the <font color = "green">fileType</font> parameter. Alternatively, if you do not have a .maf file containing annotated variants, you can set fileType = "custom" and input a custom dataframe containing columns with the names: "sample", "gene", and "variant_class".

For example: <!--Example taken from GenVisR vignette-->
```R 
# Specify the fileType 
waterfall(x=brcaMAF, fileType="MAF)

## Head the dataframes to see their contents
head(brcaMAF)
```

#### Waterfall - Identify Recurrently Mutated Genes
Due to the large number of genes mutated in at least 1 sample in this cohort, it is impossible to know which the genes actually mutated without increasing the size of the graphic device. To solve this problem, the waterfall function has 3 options. 
<ol>
<li>Use the mainRecurCutoff parameter to specify the minumum percentage of samples within the cohort that harbor a mutation in each gene. This percentage value must be represented as a decimal (e.g. 50% = 0.5).</li>
<li>The maxGenes parameter will only plot X number of the most recurrently mutated genes across the cohort.</li>
<li>The plotGenes parameter will use a character vector to plot the specified genes of interest.</li> 
</ol>

<!--Example taken from GenVisR vignette-->
```R 
# Set the mainRecurCutoff parameter
cutoff <- 0.3
waterfall(brcaMAF, mainRecurCutoff = cutoff)

# Plot the top 20 most recurrently mutated genes 
waterfall(x=brcaMAF, maxGenes = 20)

# Plot the genes TP53, RB1, BRCA1, MLL3, TTN, and MYB
genes <- c("TP53", "RB1", "BRCA1", "MLL3", "TTN", "MYB")
waterfall(x=brcaMAF, plotGenes = genes)
```

<!-- Talk about changing the order and the color -->
#### Waterfall - Customize Variant Class Order
The waterfall plot can be used to depict samples that harbor mutations in genes through other alterations, such as structural variants. These large-chromsomal aberrations include translocations (TRA), deletions (DEL), duplications (DUP), and inversions (INV). We will use structural variant data from a paper on noncirrhotic hepatocellular carcinoma to create a waterfall plot. By combining this data with the SNVs/INDELs variants, it is possible to plot the entire mutational landscape across the full sample cohort. We will use [find paper that has the sv data and the clinical data readily available].

```R
# Plot structural variant data from generic manta output
# Read in the SV data
sv_data <- 

# Assign the variant class order to include: TRA, INV, DEL, and DUP
variant_order <- c("TRA", "DEL", "DUP", "INV")
waterfall(x=sv_data, fileType="custom", variant_class_order = variant_order)

# Switch the order around to see how it changes the waterfall plot output
variant_order <- c("INV", "TRA", "DUP", "DEL")
waterfall(x=sv_data, fileType="custom", variant_class_order = variant_order)

# Change the color of the variant types 
color <- c("green", "blue", "red", "orange") ## Color order corresponds to the variant type order provided in the variable, variant_order.
waterfall(x=sv_data, fileType="custom", variant_class_order = variant_order, mainPalette = color)
```

#### Waterfall - Specify the Clinical Information
To provide both the user and reader more information on the clinical significance of the observed mutations, it is useful to illustrate relevant clinical information for each sample. This can be done by manually creating a clinical information dataframe that contains all of the variables and corresponding values for each sample. For instance, in the following example, we will want to obtain the following information: [list out the variables that we want/have]. Similar to the variant classes, both the order and the color of each clinical variables is fully customizable. 

<!--Example taken from GenVisR vignette-->
```R 
# Create clinical data
subtype <- c("lumA", "lumB", "her2", "basal", "normal")
subtype <- sample(subtype, 50, replace = TRUE)
age <- c("20-30", "31-50", "51-60", "61+")
age <- sample(age, 50, replace = TRUE)
sample <- as.character(unique(brcaMAF$Tumor_Sample_Barcode))
clinical <- as.data.frame(cbind(sample, subtype, age))

# Melt the clinical data into 'long' format.
library(reshape2)
clinical <- melt(clinical, id.vars = c("sample"))

# Run waterfall
waterfall(brcaMAF, clinDat = clinical, clinVarCol = c(lumA = "blue4", lumB = "deepskyblue", 
    her2 = "hotpink2", basal = "firebrick2", normal = "green4", `20-30` = "#ddd1e7", 
    `31-50` = "#bba3d0", `51-60` = "#9975b9", `61+` = "#7647a2"), plotGenes = c("PIK3CA", 
    "TP53", "USH2A", "MLL3", "BRCA1"), clinLegCol = 2, clinVarOrder = c("lumA", 
    "lumB", "her2", "basal", "normal", "20-30", "31-50", "51-60", "61+"))

```

#### Waterfall - Mutation Burden
The last panel of the watefall plot that has yet to be discussed is the mutational burden panel at the top. This value is determined based on the input data to the waterfall plot, and caluclates the following value in each sample: mutations per sample/coverage space * 1000000. By default, the coverage space is equivalent to the size of the SeqCap EZ Human Exome Library v2.0, but the coverage space could be changed by the user using the parameter coverageSpace. Alternatively, the user can provide a dataframe that has the claculated mutational burden for each sample. The columns for this dataframe must be "sample" and "mut_burden", and is supplied to the argument, mutBurden. 

<!--Example taken from GenVisR vignette-->
```R 
# Make the mutBurden data


```
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

## cnView
While cnSpec can illustrate the genome and cohort-wide CNV landscape, cnView focuses on copy number values for a single chromsome of a single sample. Similar to cnSpec, the values follow a relative or absolute scale specified by the <font color="green">CNscale</font> argument. The required input columns for the cnView function are chromosome, coordinate, and cn, with an optional p-value column. The chromosome of interest is specified by the <font color="green">chr</font> parameter and the <font color="green">genome</font> specification is the same as cnSpec. 

The cnView function will generate an ideogram for the identified chromsome that is on the same scale as the plot for the CNV data points. Information for the ideogram and chromosomal boundaries is obtained in the same manner as the CNspec function (from a preloaded genome or the UCSC sql database). This cytogenetic information can also be provided to the parameter <font color="green">y</font> as a dataframe with the column names: chrom, chromStart, chromEnd, name, and gieStain. <!--This overlaps significantly with the vignette--> 
The optional p-value column is provided to specify genomic coordinates in regions that are significantly amplified/deleted. If the p-value is provided in the dataframe, the transparency of the points plotted in the cnView panel will vary with low transparency occuring with lower p-values. Additionally, if you have segmentation data that contains information of the copy number loss/gain, the user can specify the segments of CNV by inputting the values as a dataframe into the <font color="green">z</font> parameter. Required columns include: chromosome, start, end, and segmean. 

 <!--Example taken from GenVisR vignette-->
```R
# create copy number data
chromosome <- "chr14"
coordinate <- sort(sample(0:106455000, size = 2000, replace = FALSE))
cn <- c(rnorm(300, mean = 3, sd = 0.2), rnorm(700, mean = 2, sd = 0.2), rnorm(1000, 
    mean = 3, sd = 0.2))
data <- as.data.frame(cbind(chromosome, coordinate, cn))

# create segment data
dataSeg <- data.frame(chromosome = c(14, 14, 14), start = coordinate[c(1, 301, 
    1001)], end = coordinate[c(300, 1000, 2000)], segmean = c(3, 2, 3))
# call cnView with included segment data
cnView(data, z = dataSeg, chr = "chr14", genome = "hg19", ideogram_txtSize = 4)
```
<!--Will make object oriented-->
## lohSpec (IN THE PROCESS OF MAKING THIS OBJECT ORIENTED)
lohSpec uses a sliding window approach to calcualte the average variant allele frequency (VAF) difference between matched tumor normal samples. Essentially, a window of a specified size (defined by <font color="green">window_size</font>) will obtain the absolute mean difference between the tumor and normal VAF (set at 0.5) for each genomic position within the window's position. This window will then move forward by, specified by the parameter <font color="green">step</font>, and again cauclate the absolute mean tumor-normal VAF difference. This value across all overlapping windows is averaged. Each value is then stored for the overlapping region and this process is repeated for each chromosome in all samples in the cohort. Plotting the X and Y chromosome is optional but can be included when a character vector containing the gender of each sample ("M" for male and "F" for female) is supplied to the <font color="green">gender</font> parameter. 

<!--Example taken from GenVisR vignette-->
```R
## Generate LOHspec with default parameters
lohSpec(x=HCC1395_Germline)
````

<!--Talk about caveats with LOH-->
When interpreting this visualization and other LOH data, it is important to understand that not all LOH corresponds to an actual deletion event. For instance, an amplification of a specific allele/genomic region can appear as LOH since one region will have significantly more copies than the corresponding region on the other chromosome. Additionally, copy neutral loss of heterozygosity can occur where an individual receives a chromosomal pair from one parent and nothing from the other parent. While this situation may appear as LOH since there is only a single chromosomal type present, the individual still has 2 chromosomes (hence neutral copy number). Therefore, it is important to look at coverage to see LOH truly results in a deleted region. 

## genCov


## More GenVisR Functions
For a full list of GenVisR functions with descriptions and examples, see the GenVisR reference manual and vignette on the Bioconductor page at [https://bioconductor.org/packages/release/bioc/html/GenVisR.html](https://bioconductor.org/packages/release/bioc/html/GenVisR.html). 
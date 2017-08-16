---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to waterfall plots
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-02-01
---

After sequencing a set of samples a common question to ask is what mutations are present in my cohort? This type of data is often displayed in a heatmap like structure with rows and columns denotating genes and samples with these variables ordered in a hierarchical pattern based on recurrence. Such plots make viewing mutually exclusive and co-occurring mutational event a trivial task. The [waterfall()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/waterfall) function from the [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) package makes it easy to create these types of charts while also allowing additional data in the context of sample, and gene data to be added to the plot. In this tutorial we will use the [waterfall()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/waterfall) function to re-create panel C of figure 2 from the paper ["A Phase I Trial of BKM120 (Buparlisib) in Combination with Fulvestrant in Postmenopausal Women with Estrogen Receptor-Positive Metastatic Breast Cancer."](https://www.ncbi.nlm.nih.gov/pubmed/26563128).

{% include figure.html image="/assets/GenVisR/BKM120_Waterfall_Final.png" width="950" %}

### Installing GenVisR and loading Data
To begin let's install and load the [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html) library from bioconductor. We will also need to load the mutation data from [Supplemental Table 3](http://clincancerres.aacrjournals.org/content/suppl/2015/11/12/1078-0432.CCR-15-1745.DC1), a version that has been converted to a tab delimited format is available for download [here](http://genomedata.org/gen-viz-workshop/GenVisR/BKM120_Mutation_Data.tsv). We also need some supplemental information regarding the clinical data which is available [here](http://genomedata.org/gen-viz-workshop/GenVisR/BKM120_Clinical.tsv) and mutation burden data available [here](http://genomedata.org/gen-viz-workshop/GenVisR/BKM120_MutationBurden.tsv).

```R
# install and load GenVisR
source("https://bioconductor.org/biocLite.R")
biocLite("GenVisR")
library(GenVisR)

# load relevant data from the manuscript
mutationData <- read.delim("BKM120_Mutation_Data.tsv")
clinicalData <- read.delim("BKM120_Clinical.tsv")
mutationBurden <- read.delim("BKM120_MutationBurden.tsv")
```

# creating the inital plot
The [waterfall()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/waterfall) function is designed to work with specific filetypes read in as data frames, the default being [MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) files, however the option exist to use custom file types as long as the column names "sample", "gene", and "variant_class" are present. This is done by setting `fileType="Custom"`. Let's go ahead and re-format our mutation data and create an initial plot using this parameter. Note that when using a custom file type the priority of mutation types must be specified with the `variant_class_order` parameter which accepts a character vector. We'll explain what this does a bit more in the next section, for now just give it a character vector containing all mutation types present in `mutationData`.
```R
# reformat the mutation data for waterfall()
mutationData <- mutationData[,c("patient", "gene.name", "trv.type")]
colnames(mutationData) <- c("sample", "gene", "variant_class")

# create an inital plot
waterfall(mutationData, fileType = "Custom", variant_class_order=as.character(unique(mutationData$variant_class)))
```

# Setting the mutation hierarchy
As you can see from the previous plot, [waterfall()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/waterfall) plots mutations in a method similar to a heatmap, however the astute user may be asking what happens in cases where there are two different mutations types for the same cell in the heatmap. The answer is [waterfall()](https://www.rdocumentation.org/packages/GenVisR/versions/1.0.4/topics/waterfall) picks one based on a mutation hierarchy. For a known file format such as [MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) this hierarchy is pre set such that the more deleterious mutations are higher on the hierarchy however for a "custom" file type there is no apriori knowledge of what this order should be and thus must be set by the user. As mentioned this is done with the `variant_class_order` parameter which accepts a character vector of all mutation types within the mutation data input in order of most to least deleterious. The effect this parameter has is best exemplified in the example below.
```R
# make sure seed is set to 426 to reproduce!
set.seed(426)

# install and load the gridExtra package
#install.packages("gridExtra")
library(gridExtra)

# Create a data frame of random elements to plot
inputData <- data.frame(sample = sample(letters[1:5], 20, replace = TRUE), gene = sample(letters[1:5], 20, replace = TRUE), variant_class = sample(c("x", "y", "z"), 20, replace = TRUE))

# choose the most deleterious to plot with y being defined as the most
# deleterious
most_deleterious <- c("y", "z", "x")

# plot the data with waterfall using the 'Custom' parameter
p1 <- waterfall(inputData, fileType = "Custom", variant_class_order = most_deleterious, mainXlabel = TRUE, out="grob")

# change the most deleterious order
p2 <- waterfall(inputData, fileType = "Custom", variant_class_order = rev(most_deleterious), mainXlabel = TRUE, out="grob")

# arrange the two plots side by side
grid.arrange(p1, p2, ncol=2)
```

{% include figure.html image="/assets/GenVisR/waterfall_hierarchy_example.png" width="950" %}

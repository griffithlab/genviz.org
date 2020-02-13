---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Data Munging with Data.table
categories:
    - Module-02-R
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-02
---

In the previous section we went through the basics of R, that will get you pretty far but often you will need to manipulate your data in some way to answer a biological question. In the R ecosystem there are three well known competing libraries in order to help in this regard. The first are the base R functions which we touched on a bit previously such as `[]` to manipulate a data.frame and aggregate() in order to apply some function to the data. A second option is the [dplyr](https://dplyr.tidyverse.org/), part of the tidyverse which is much faster and the base R functions and has what many find to be a more intuitive syntax. The third option is the [data.table](https://rdatatable.gitlab.io/data.table/), this is what we will be going over in this section. It is extremly memory effecient and is consistently among the fastest solutions.

Let's start by loading the library into R and reading in some data. To do this we will use fread() instead of the normal read.delim() fuction.

```R
# load package
install.packagees("data.table")
library(data.table)

# read data
varData <- fread("http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv")
class(varData)
```

If you were paying attention when you loaded the data.table library you might have noticed that data.table automatically recognized that your computer had multiple cores and registered them automatically. So when possible data.table will use multi-threading without you having to to anything. You'll also note that the class of varData is both a data.table and data.frame. A nice thing about the data.table package is that all functions that work with data.frames should also work with the data.table class as well due to inheritance.

Now that we have are data read in let's go over how to subset rows from the data.table structure. This will work mostly like how we would subset data.frames with a couple of key differences. 1. data.table will automatically assume you want to subset rows if there is no comma in the square brackets and 2. you can just supply the column name to data.table to subset via a boolean expression. Let's have a look with some examples

```R
# subset by row index
varData[1:5]

# extract all columns where chromosome is 22
varData[chromosome_name == 22]
```

We can also select columns via a similar method to how we're used to in data.frame, using character vectors or column indices, however the best practice is to simply use a list. the .() function in data.table is an alias for list().

```R
# select out the chromosome_name column
varData[,2]
varData[,"chromosome_name"]
varData[,.(chromosome_name)]
```

In addition to selecting column names we can compute results directly inside a data.table from within the .() list of column names.

```R
# find the mean tumor vaf
varData[,.(mean(tumor_VAF))]

# equivalent to above but we put the result in a named column
varData[,.(tumor_vaf = mean(tumor_VAF))]

# find the mean tumor vaf and mean tumor coverage
varData[,.(tumor_mean_vaf = mean(tumor_VAF), tumor_mean_cov=mean(tumor_var_count + tumor_ref_count))]

# we can apply this to select rows as well
varData[sample=="H_ML-08-0075-001-1127127",.(tumor_vaf = mean(tumor_VAF))]
```

At this time it is appropriate to introduce a special symbol within the data.table package, the .N. It stands for the total number of rows within the data.table, essentially equivalent to the nrow() function. We can use this function to count filtered results.

```R
# find how many entries chromosome 1 has
varData[chromosome_name=="1",.N, by=.(chromosome_name)]
```

Thus far we have dealt with topics that have some simalarity with data.frames, however let's now dive into something unique to the data.table structure, grouping within data.tables. Essentially we can perform computations on columns by groups of values. For those familiar with base R this replaces the aggregate() function and can be used within a data.table call. Let's look at some examples!

```R
# count the number of entries per chromosome
varData[,.N,by=.(chromosome_name)]

# count the number of entries for strand and chromosome
varData[,.N,by=.(chromosome_name, strand)]

# count the number of entries per chromosome for the discovery set
varData[dataset=="discovery",.N,by=.(chromosome_name)]
```

We've seen some examples on manipulating data.table objects so at this point let's talk about a best practice feature of data.table, the ability to update a data.table object by reference. In current versions of R, 3.X.X at the time of this writing, when manipulating a column of a data frame the entire column get's copied into internal memory, the column is updated, and is then added into the original data frame. This means that if you had a single column data frame taking up 2GB of memory and wanted to update it that operation would use 4GB. This is part of the reason R has a reputation as a memory hog. Fortunatley data.table has the option of updating columns by reference with the `:=` operator. When using this no copy of the data is made but rather a reference is made mapping the old value to the new value. Let's go over some examples for how it works. We will used the gc() function for a rough estimate of memory ussage as the ussual methods don't account for garbage collection which would skew results in this case. See this [link](https://stackoverflow.com/questions/58278838/memory-profiling-with-data-table) for an explanation as to why.

```R
# create a data.table and data.frame for illustration purposes
myDF <- data.frame("V1"=c(1:10000), "V2"=c(1:10000))
myDT <- data.table("V1"=c(1:10000), "V2"=c(1:10000))

# profile the data.frame and data.table memory ussage for adding two columns
gc(myDF$V3 <- myDF$V1 + myDF$V2)
gc(myDT[,"V3" := V1 + V2])

# modify just a single value and profile
gc(myDF[1,"V1"] <- 100)
gc(myDT[1,"V1" := 100])

# did I mention you can modify by reference multiple columns at once
myDT[,c("rev_v1", "rev_v2") := .(rev(V1), rev(V2))]
```

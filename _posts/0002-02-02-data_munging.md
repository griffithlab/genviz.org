---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Data Munging with Data.table
categories:
    - Module-02-R
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-02
---

In the previous section we went through the basics of R, that will get you pretty far but often you will need to manipulate your data in some way to answer a biological question. In the R ecosystem there are three well known competing libraries in order to help in this regard. The first are the base R functions which we touched on a bit previously such as `[]` to manipulate a data.frame and aggregate() in order to apply some function to the data. A second option is the [dplyr](https://dplyr.tidyverse.org/), part of the tidyverse which is much faster than the base R functions and has what many find to be a more intuitive syntax. The third option is the [data.table](https://rdatatable.gitlab.io/data.table/) library, this is what we will be going over in this section. It is extremly memory effecient and is consistently among the fastest solutions. throughout this course the base-r way will be shown alongside the data.table for comparison of syntax, these lines will be denoted with `#!`.

Let's start by loading the library into R and reading in some data. To do this we will use fread() instead of the normal read.delim() fuction.

```R
# load package
install.packages("data.table")
library(data.table)

# read data
#! varDataDF <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv")
varDataDT <- fread("http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv")
class(varData)
```

If you were paying attention when you loaded the data.table library you might have noticed that data.table automatically recognized that your computer had multiple cores and registered them automatically. So when possible data.table will use multi-threading without you having to to anything. You'll also note that the class of varData is both a data.table and data.frame. A nice thing about the data.table package is that all functions that work with data.frames should also work with the data.table class as well due to inheritance.

Before we really get started we should go over the overarching theme of data.table syntax, `DT[i,j,by]`. The `i` argument similar to dataframes, are on which rows to act. The `j` argument are what to do on columns, unlike dataframes we can do more than just select columns. Finally the `by` argument stands for group by what, something unique to data.tables. As we move along this syntax should start to make more sense.

Now let's go over how to subset rows from the data.table structure. This will work mostly like how we would subset data.frames with a couple of key differences. 1. data.table will automatically assume you want to subset rows if there is no comma in the square brackets and 2. you can just supply the column name to data.table to subset via a boolean expression. Let's have a look with some examples

```R
# subset by row index
#! varDataDF[1:5,]
varDataDT[1:5]

# extract all columns where chromosome is 22
#! varDataDF[varDataDF$chromosome_name == 22,]
varDataDT[chromosome_name == 22]
```

We can also select columns via a similar method to how we're used to in data.frame, using character vectors or column indices, however the best practice is to simply use a list. the .() function in data.table is an alias for list().

```R
# select out the chromosome_name column
#! varDataDF[,2]
varDataDT[,2]
#! varDataDF[,"chromosome_name"]
varDataDT[,"chromosome_name"]
#! varDataDT[,"chromosome_name"]
varDataDT[,.(chromosome_name)]
```

In addition to selecting column names we can compute results directly inside a data.table from within the .() list of column names.

```R
# find the mean tumor vaf
#! mean(varDataDF$tumor_VAF)
varDataDT[,.(mean(tumor_VAF))]

# equivalent to above but we put the result in a named column
#! data.frame("tumor_vaf"=mean(varDataDF$tumor_VAF))
varDataDT[,.(tumor_vaf = mean(tumor_VAF))]

# find the mean tumor vaf and mean tumor coverage
#! data.frame("tumor_mean_vaf"=mean(varDataDF$tumor_VAF), "tumor_mean_cov"=mean(varDataDF$tumor_var_count + varDataDF$tumor_ref_count))
varDataDT[,.(tumor_mean_vaf = mean(tumor_VAF), tumor_mean_cov=mean(tumor_var_count + tumor_ref_count))]

# we can apply this to select rows as well
#! data.frame("tumor_vaf"=mean(varDataDF[varDataDF$sample=="H_ML-08-0075-001-1127127","tumor_VAF"]))
varDataDT[sample=="H_ML-08-0075-001-1127127",.(tumor_vaf = mean(tumor_VAF))]
```

At this time it is appropriate to introduce a special symbol within the data.table package, the .N. It stands for the total number of rows within the data.table, essentially equivalent to the nrow() function. We can use this function to count filtered results.

```R
# find how many entries chromosome 1 has
#! table(varDataDF[varDataDF$chromosome_name == "1","chromosome_name"], exclude = NULL)
varDataDT[chromosome_name=="1",.N, by=.(chromosome_name)]
```

Thus far we have dealt with topics that have some simalarity with data.frames, however let's now dive into something unique to the data.table structure, grouping within data.tables. Essentially we can perform computations on columns by groups of values. For those familiar with base R this replaces the aggregate() function and can be used within a data.table call. Let's look at some examples!

```R
# count the number of entries per chromosome
#! as.data.frame(table(varDataDF[,"chromosome_name"], exclude = NULL))
varDataDT[,.N,by=.(chromosome_name)]

# count the number of entries for strand and chromosome
#! as.data.frame(table(varDataDF[,c("chromosome_name", "strand")], exclude = NULL))
varDataDT[,.N,by=.(chromosome_name, strand)]

# count the number of entries per chromosome for the discovery set
#! as.data.frame(table(varDataDF[varDataDF$dataset=="discovery",c("chromosome_name")], exclude = NULL))
varDataDT[dataset=="discovery",.N,by=.(chromosome_name)]
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

Now that we've gone over the very basics for data.table i.e. the `DT[i,j,by]` syntax and how to assign columns on the fly, let's introduce some more concepts to get a feel for exactly how powerfull data.table can be. To begin, let's introduce the other special variables in data.table starting with `.GRP`. `.GRP` holds the grouping id from the `by` argument. Let's say for for argument a bug was introduced in the code causing the variants on even chromosomes for the extension cohort to be 1 base off. How could you fix this in base-r succinctly, I actual don't have a succinct answer for base R, but with data.table it's fairly straight forward.

```R
# see below for an explanation
varDataDT[order(chromosome_name, dataset), "group":=.GRP, by=.(chromosome_name, dataset)][group %% 4 == 0,"new_start":=start + 1]
```

So this is unlikely to happen but does allow us to introduce a couple concepts, some of which are new. To start things off we first act on rows by ordering by "chromosome_name" and "dataset" via `order(chromosome_name, dataset)` we also tell data.table to group by "chromosome_name" and "dataset" `by=.(chromosome_name, dataset)`, and to assign that grouping to a new column called "group" `"group":=.GRP`. At this point every group which is divisible by 4 contains the value we wish to adjust so we chain the data.table we jsut created to another data.table expression simply by adding `[]` brackets. From there we can simply filter to groups which are divisible by 4 using the modulus operator `group %% 4 == 0` and make a new column increasing the start position by 1 `"new_start":=start + 1`.

Okay, I know unrealistic example, but you get the idea, try doing what we did above with base-r code. Let's go over another special variable with a more realistic application. The `.SD` variable stands for subset of data and essentially stores the subsets of data from the `by` argument as data.tables. It is commonly used with thhe `.SDcols` variable which specifys the columns to return for the data subsets. Let's imagine that for each sample you need to return the first variant in varDataDT. Here's the data.table way to do it.

```R
varDataDT[order(chromosome_name, start), .SD[1], by=.(Simple_name), .SDcols=c("Simple_name", "chromosome_name", "start")]
```

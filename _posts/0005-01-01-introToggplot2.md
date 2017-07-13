---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to ggplot2
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-01-01
---

There are three primary graphics program available within the R environment. Base-r graphics are installed by default and provide a simple mechanism to quickly create graphs, lattice is a graphics program influenced by trellis graphics, and ggplot2 is a graphics program based on the grammar of graphics idealogy. In this course we will be focusing on ggplot2 as our graphics package of choice and use it to explore Supplemental Table S5 of the paper ("Recurrent somatic mutations affecting B-cell receptor signaling pathway genes in follicular lymphoma")[http://www.bloodjournal.org/content/129/4/473/tab-figures-only].

### basic ggplot2 syntax
ggplot is based on a system of layering graphical object to create a final plot and takes data frames as it's input. Let's start with loading the ggplot2 library and plotting the tumor variant allele fractions and coverage.

```R
# install the ggplot2 library and load it
install.packages("ggplot2")
library(ggplot2)

# load Supplemental Table S5
variantData <- read.delim("~/Desktop/ggplot2ExampleData.tsv.txt")

# make a coverage column since this doesn't exist yet
variantData$tumor_COV <- variantData$tumor_ref_count + variantData$tumor_var_count

# start making the plot
p1 <- ggplot(data=variantData, aes(x=tumor_VAF, y=tumor_COV))
p1
```

### geom and aes
you should see a blank plot at this point so what is going on. Well we invoked ggplot with ggplot() and gave it the data frame we wish to plot. We then supplied aesthetic mappings with aes() which told ggplot which columns should be assigned to the geometric objects aesthetics. In this specific case we are telling ggplot that the data is in "varaintData" and we want the column tumor_VAF to be assigned to the x-axis and tumor_COV to be assigned to the y-axis. ggplot has done this and we can see it has automatically determined the axis scales given the ranges in the data supplied however nothing is actually plotted because we have not said what geometric object to use. let's try adding a geometric object as a layer to the plot and see what happens, all layers can be added with a "+".

```R
# add a point geom to the plot (method 1)
p1 <- p1 + geom_point()
p1

# this is equivalent to above (method 2)
p2 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=tumor_COV))
p2

# illustrating the difference between continous and discrete data
p3 <- ggplot(data=variantData, aes(x=chromosome_name, y=tumor_VAF), position="jitter")
```

At this point you should see a plot with points showing the variant allele fraction distribution in relation to the corresponding read coverage. While method 1 and method 2 produce identical plots we should briefly talk about the differences between the two. In method 1 we invoke ggplot with a data frame and aesthetic mapping, this information is carried over to all subsequent layers where appropriate. Conversley in method 2 we are invoking ggplot but defining the data and mapping only within the point geometric object, this is useful when plotting different datasets on the same plot as we will see later. We should also note that geoms can behave differently depending on if the data supplied to them is continuous or discrete. In the example above (plot p3), we can see that the points have been binned  by chromosome name on the x-axis while the contiuous scale is kept on the y-axis. we also told ggplot to jitter the points along the x-axis to provide better resolution.

### axis scaling and manipuation
Going back to our original example (p1), the majority of points look like they have a coverage < 500 however the coverage outliers in this data set is causing that area which is perhaps the most interesting to take up a relatively small proportions of the plot. There are a few ways we can provide more resolution in this area, lets' test a few out and discuss them.

```R
# method 1, set y limits for the plot
p4 <- p1 + scale_y_continous(limits=c(0, 500))
p4 <- p1 + ylim(c(0, 500)) # shortcut method
p4

# method 2, transform the actual values in the data frame
p5 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=log2(tumor_COV)))
p5

# method 3, transform the axis using ggplot
p6 <- p1 + scale_y_continuous(trans="log2")
p6
```

Perhaps the easiest thing to do is simply only plot the points within the range we want, we can achieve this by adjusting the [scale_y_continuous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html) layer. [ylim()](http://ggplot2.tidyverse.org/reference/lims.html) is a shortcut which doesn the same thing. You'll see a warning when doing this, all it's saying is that it's removing rows from the data frame not within the specified range. There is an "out of bounds" parameter within scale_y_continuous() to control what happens with these points but perhaps this isn't the best method. In method 2 we actually transform the data values themselves altering applying a log2 transform. While this works beautifully it perhaps isn't intuitive to interpret the log2 of coverage. A better method might be to not transform the values but rather adjust the plotting scale itself (method 3).

This plot is looking pretty good but is lacking color we can add some by simply defining that aesthetic.

```R
# add some color to the plot
p1 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=tumor_COV, colour="red")) + scale_y_continuous(trans="log2")

p1 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=tumor_COV), colour="red") + scale_y_continuous(trans="log2")

```

### displaying additional variables

### faceting

### wide vs long format

### ggvis
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### additional resources
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

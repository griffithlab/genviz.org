---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to ggplot2
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-01-01
---

There are three primary graphics program available within the R environment. base-r graphics are installed by default and provide a simple mechanism to quickly create graphs. lattice is a graphics program influenced by trellis graphics. ggplot2 is a graphics program based on the grammar of graphics idealogy. In this course, we will be focusing on ggplot2 as our graphics package of choice and use it to explore Supplemental Table S5 of the paper ["Recurrent somatic mutations affecting B-cell receptor signaling pathway genes in follicular lymphoma"](http://www.bloodjournal.org/content/129/4/473/tab-figures-only).

### Basic ggplot2 syntax
ggplot is based on a system of layering graphical objects to create a final plot, utilizing data frames as its input. We will start by installing and loading the ggplot2 library. After importing our data ('ggplot2ExampleData.tsv.txt'), we will modify this data frame to include a 'coverage' (tumor_COV) variable. Then we can call the variantData data frame in our ggplot() function and compare the coverage variable to the variant allele frequency (tumor_VAF).

```R
# install the ggplot2 library and load it
install.packages("ggplot2")
library(ggplot2)

# load Supplemental Table S5
variantData <- read.delim("ggplot2ExampleData.tsv")

# make a coverage column since this doesn't exist yet
variantData$tumor_COV <- variantData$tumor_ref_count + variantData$tumor_var_count

# start making the plot
p1 <- ggplot(data=variantData, aes(x=tumor_VAF, y=tumor_COV))
p1
```

### geom and aes
The variable p1 generates a blank plot in the bottom right "Plot" window in Rstudio. We invoked ggplot with [ggplot()](http://ggplot2.tidyverse.org/reference/ggplot.html) and specified the data frame we are trying to plot. We then supplied aesthetic mappings with [aes()](http://ggplot2.tidyverse.org/reference/aes.html), which told ggplot which columns should be assigned to the geometric objects aesthetics. In this specific case, we are telling ggplot that the data is in the data frame "variantData", the column tumor_VAF should be plotted along the x-axis, and tumor_COV should be plotted along the y-axis. ggplot has determined the axis scales, given the ranges in the data supplied. However, nothing is actually plotted because we have not specified what geometric object to use. The geometric objects in ggplot describe how we are plotting the data (e.g. geom_point() for scatterplots, geom_bar() for bar charts). Geometric objects are added as plot layers to the ggplot() command using a [+](http://ggplot2.tidyverse.org/reference/gg-add.html).

```R
# add a point geom to the plot (method 1)
p1 <- p1 + geom_point()
p1

# this is equivalent to above (method 2)
p2 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=tumor_COV))
p2

# illustrating the difference between continous and discrete data
p3 <- ggplot() + geom_point(data=variantData, aes(x=chromosome_name, y=tumor_VAF), position="jitter")
p3
```

At this point you should see a plot with points showing the variant allele fraction distribution in relation to the corresponding read coverage. While method 1 and method 2 produce identical plots we should briefly talk about the differences between the two. In method 1 we invoke ggplot with a data frame and aesthetic mapping, this information is carried over to all subsequent layers where appropriate. Conversley in method 2 we are invoking ggplot but defining the data and mapping only within that geometric object, this is useful when plotting different datasets on the same plot as we will see later. We should also note that geoms can behave differently depending on if the data supplied to them is continuous or discrete. In the example above (plot p3), we can see that the points have been binned  by chromosome name on the x-axis while the contiuous scale is kept on the y-axis. we also told ggplot to jitter the points along the x-axis to provide better resolution. A complete list of available geoms within ggplot is available [here](http://ggplot2.tidyverse.org/reference/#section-layer-geoms).

### axis scaling and manipuation
Going back to our original example (p1), the majority of points look like they have a coverage < 500 however there are outliers in this data set which is causing the area which is perhaps the most interesting to take up a relatively small proportions of the plot. There are a few ways we can provide more resolution in this area, lets' test a few out and discuss them.

```R
# method 1, set y limits for the plot
p2 <- p1 + scale_y_continuous(limits=c(0, 500))
p2 <- p1 + ylim(c(0, 500)) # shortcut method
p2

# method 2, transform the actual values in the data frame
p3 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=log2(tumor_COV)))
p3

# method 3, transform the axis using ggplot
p1 <- p1 + scale_y_continuous(trans="log2")
p1
```

Perhaps the easiest thing to do is simply plot just the points within the range we want, we can achieve this by adjusting the [scale_y_continuous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html) layer. [ylim()](http://ggplot2.tidyverse.org/reference/lims.html) is a shortcut which does the same thing. You'll see a warning when doing this, all it's saying is that it's removing rows from the data frame not within the specified range. There is an "out of bounds" parameter within [scale_y_continuous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html) to control what happens with these points but perhaps this isn't the best method to achieve our goal. In method 2 we actually transform the data values themselves applying a log2 transform. While this works beautifully it perhaps isn't intuitive to interpret the log2 of read coverage. A better method might be to not transform the values but rather adjust the plotting scale itself (method 3).

### using different aesthetics

This plot is looking pretty good but is lacking color we can add some by simply defining that aesthetic. We can specify a color by either the hex code or by naming it from R's internal colour pallette, a full list of which is available [here](http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf), alternatively you can list colors by typing colors() in the R terminal.

```R
# list colors in R
colors()

# what happens when we try and add colour withing the aesthetic
p2 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=tumor_COV, colour="#68228B")) + scale_y_continuous(trans="log2")
p2

# and what happens when we try and add colour within the geom
p3 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=tumor_COV), colour="#68228B") + scale_y_continuous(trans="log2")
p3

```

Above we chose "darkorchid4" which has a hex value of "#68228B", however the points in the first plot (p2) are red not purple as would be expected so what is going on? The answer is reliant on how ggplot interprets aesthetic mappings. When we specified our quoted hex code in the aesthetic ggplot assummed we wanted to add another variable to the data. It did this for us and then used it's internal color scheme to color that variable. To actually color the points purple we can specify the color outside the aesthetic mapping as in the latter plot (p3).

We can use this to answer a few questions, what if we wanted to know if the "discovery" or "extension" cohorts within the data have a higher tumor purity. Let's use [geom_density()](http://ggplot2.tidyverse.org/reference/geom_density.html), which will plot a density kernel, in combination with color to find out.

```
# get a density curve of tumor vafs
p1 <- ggplot() + geom_density(aes(x=tumor_VAF, colour=dataset))
p1

# let's add a bit more detail
p2 <- ggplot() + geom_density(data=variantData, aes(x=tumor_VAF, fill=dataset), alpha=.75, colour="black", adjust=.5)
p2

# and let's change the colors some more
p3 <- p2 + scale_fill_manual(values=c("discovery"="#a13242", "extension"="#1a2930"))
p3
```

In the first plot we've told ggplot to differentiate the data based on the "dataset" column using color. As we can see it separated the data out based on this variable resulting in two density curves. In the second plot we are telling ggplot to do the same thing however we are using the aesthetic fill to differentiate the data instead and are globally assigning the colour and alpha which correspond to the transparency and line colour respectively. We also set the adjust parameter which will reduce the smoothing [geom_density()](http://ggplot2.tidyverse.org/reference/geom_density.html) uses when computing it's estimate. Finally in the last plot we manually adjust the fill colours ggplot uses with [scale_fill_manual()](http://ggplot2.tidyverse.org/reference/scale_manual.html)

### faceting

Depending on the geometric object used there are up to 10 ways to map an aesthetic to a variable. These are with the x-axis, y-axis, fill, colour, shape, alpha, size, labels, and facets. Faceting in ggplot allows us to quickly create multiple related plots at once with a single command. Let's try and answer a few quick questions about our data using facets.

### displaying additional variables

### wide vs long format

### ggvis
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### additional resources
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

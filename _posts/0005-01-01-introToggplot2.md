---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to ggplot2
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-01-01
---

There are at least three primary graphics programs available within the R environment. A package for [base R graphics](https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/00Index.html) is installed by default and provides a simple mechanism to quickly create graphs. [lattice](https://cran.r-project.org/web/packages/lattice/index.html) is another graphics package that attempts to improve on base R graphics by providing better defaults and the ability to easily display multivariate relationships. In particular, the package supports the creation of trellis graphs - graphs that display a variable or the relationship between variables, conditioned on one or more other variables. Finally, [ggplot2](http://ggplot2.org/) is a graphics program based on the grammar of graphics idealogy, and will be the primary focus of this course. 

In this module, we will explore basic use of ggplot2 to plot genomic data. For illustration, we will use a set of mutation data from Supplemental Table S5 of the paper ["Recurrent somatic mutations affecting B-cell receptor signaling pathway genes in follicular lymphoma"](http://www.bloodjournal.org/content/129/4/473/tab-figures-only). You can download a cleaned up version of Supplemental Table S5 at [http://www.genomedata.org/gen-viz-workshop/ggplot2ExampleData.tsv](http://www.genomedata.org/gen-viz-workshop/ggplot2ExampleData.tsv)

### Wide vs long format
Before we begin it is important to know that ggplot expects the data passed to it to be of class data.frame. Further the data should be in long instead of wide format. This simply means that instead of each non-id variable having it's own column there should be a column/columns designating a key/value pair. We can change between wide and long formats with the [dcast()](https://www.rdocumentation.org/packages/reshape2/versions/1.4.2/topics/cast) and [melt()](https://www.rdocumentation.org/packages/reshape2/versions/1.4.2/topics/melt) functions from the reshape2 package. For simplicity our example dataset [ggplot2ExampleData.tsv](http://www.genomedata.org/gen-viz-workshop/ggplot2ExampleData.tsv) is already in long format.

{% include figure.html image="/assets/long_v_wide.png" width="750" %}

### Introducing ggplot2 syntax
ggplot is based on a system of layering graphical objects to create a final plot, and as mentioned utilizes data frames as its input. We will start by installing and loading the [ggplot2](http://ggplot2.tidyverse.org/) library. After importing our data ('ggplot2ExampleData.tsv'), we will modify this data frame to include a 'coverage' (tumor_COV) variable. Then we can call the variantData data frame in our [ggplot()](http://ggplot2.tidyverse.org/reference/ggplot.html) function and compare the coverage variable to the variant allele frequency (tumor_VAF).

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

### Geometric objects and aesthetic mapping
The object stored in variable p1 will generate a blank plot in the bottom right "Plots" window of Rstudio. We invoked ggplot with the function [ggplot()](http://ggplot2.tidyverse.org/reference/ggplot.html) and specified the data frame we are trying to plot. We then supplied aesthetic mappings with [aes()](http://ggplot2.tidyverse.org/reference/aes.html). In essence this is specifying which columns ggplot should assign to the geometric objects aesthetics. In this specific case, we are telling ggplot that the data is in the data frame "variantData", the column tumor_VAF should be plotted along the x-axis, and tumor_COV should be plotted along the y-axis. ggplot has determined the axis scales, given the ranges in the data supplied. However you will notice that nothing is actually plotted, this is because we have not specified what [geometric object](http://ggplot2.tidyverse.org/reference/index.html#section-layer-geoms) to use. At this point we have passed what data we want plotted however we have not specified how it should be plotted. This is what the various [geometric objects](http://ggplot2.tidyverse.org/reference/index.html#section-layer-geoms) in ggplot are used for (e.g. [geom_point()](http://ggplot2.tidyverse.org/reference/geom_point.html) for scatterplots, [geom_bar()](http://ggplot2.tidyverse.org/reference/geom_bar.html) for bar charts). To do this [geometric object](http://ggplot2.tidyverse.org/reference/index.html#section-layer-geoms) are added as plot layers to the [ggplot()](http://ggplot2.tidyverse.org/reference/ggplot.html) base command using a [+](http://ggplot2.tidyverse.org/reference/gg-add.html).

```R
# add a point geom to the plot (method 1)
p1 <- p1 + geom_point()
p1

# this is equivalent to above (method 2)
p2 <- ggplot() + geom_point(data=variantData, aes(x=tumor_VAF, y=tumor_COV))
p2
```

Both plot p1 and plot p2 generate a scatter plot, comparing the column tumor_COV on the y-axis to the column tumor_VAF on the x-axis. While the plots generated by the p1 and p2 variables may appear identical, we should briefly address their differences. In method 1 (plot p1), we invoke [ggplot()](http://ggplot2.tidyverse.org/reference/ggplot.html) with a data frame (variantData) and an aesthetic mapping of tumor_VAF and tumor_COV for x and y respectively. In this method the information is passed to all subsequent [geometric objects](http://ggplot2.tidyverse.org/reference/index.html#section-layer-geoms) and is used as appropriate in those object. In this case, the only geometric object we include is [geom_point()](http://ggplot2.tidyverse.org/reference/geom_point.html). The [geom_point()](http://ggplot2.tidyverse.org/reference/geom_point.html) layer is then added using the information passed from [ggplot()](http://ggplot2.tidyverse.org/reference/ggplot.html). Conversely, in method 2 (plot p2), we invoke [ggplot()](http://ggplot2.tidyverse.org/reference/ggplot.html) without defining the data or aesthetic mapping. This information is specified directly in the [geom_point()](http://ggplot2.tidyverse.org/reference/geom_point.html) layer. If any additional geometric objects were added as layers to the plot, we would specifically have to define the data and aesthetics within each additional layer. This is especially useful when plotting multpile datasets on the same plot (we will explore this later on).

We should also note that geometric objects can behave differently, depending upon whether the plotted variables are continuous or discrete. In the example below (plot p3), we can see that the points have been binned by chromosome name on the x-axis, while the numeric values sorted in the column "tumor_VAF" are plotted along a continuous y-axis scale. The position of the points (specified by position='jitter' in the [geom_point()](http://ggplot2.tidyverse.org/reference/geom_point.html) object) shifts the points apart horizontally to provide better resolution. A complete list of available geoms within ggplot is available [here](http://ggplot2.tidyverse.org/reference/#section-layer-geoms).

```R
# illustrating the difference between continous and discrete data
p3 <- ggplot() + geom_point(data=variantData, aes(x=chromosome_name, y=tumor_VAF), position="jitter")
p3
```

### Axis scaling and manipulation
Going back to our original example (plot p1), the majority of points look like they have a coverage < 500x. However, there are outliers in this data set causing the majority of points to take up a relatively small portion of the plot. We can provide more resolution to this by the following methods:
1. limiting the y-axis scale using [scale_y_continous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html) or [ylim()](http://ggplot2.tidyverse.org/reference/lims.html)
2. transforming the numeric values by [log2()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/log) on the y-axis.
3. transforming the y-axis to a log2 scale by specifying trans within the [scale_y_continous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html) layer.

Note that these transformations can be applied to the x axis as well ([scale_x_continous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html), [xlim()](http://ggplot2.tidyverse.org/reference/lims.html), etc.), as long as the x-axis is mapped to data appropriate for a continuous scale.

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

Note that adjusting the [scale_y_continuous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html) layer will plot only the points within the specified range by setting the limits. [ylim()](http://ggplot2.tidyverse.org/reference/lims.html) is a shortcut that achieves the same thing. You'll see a warning when doing this, stating that rows are being removed from the data frame that contain values outside of the specified range. There is an "out of bounds" parameter within [scale_y_continuous()](http://ggplot2.tidyverse.org/reference/scale_continuous.html) to control what happens with these points, but this isn't necessarily the best method for this particular plot. In method 2 (plot p2), we actually transform the data values themselves by applying a log2 transform. This method allows us to better visualize all of the points, but it is not intuitive to interpret the log2 of a value (tumor coverage). Alternatively, method 3 (plot p3) does not transform the values, but adjusts the scale the points are plotted on and does a better job of displaying all of our data without having to convert the log2 values.

### Applying different aesthetics
While these plots look pretty good, we can make them more aesthetically pleasing by defining the color of the points within the aesthetic. We can specify a color by either the hex code or by naming it from R's internal color pallette, a full list of which is available [here](http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf). Alternatively, you can list colors by typing [colors()](https://www.rdocumentation.org/packages/grDevices/versions/3.4.1/topics/colors) in the R terminal.

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

Above we chose "darkorchid4" which has a hex value of "#68228B". However the points in the first plot (p2) are red and not the expected purple color. Our points are appearing miscolored based upon how ggplot is interpreting the aesthetic mappings. When we specified our quoted hex code in the aesthetic in p2, ggplot assummed we wanted to add another variable to the data. It did this for us and then used its internal color scheme to color that variable. By specifying the color outside the aesthetic mapping, geom_point knows to apply the color 'darkorchid4' to all of the points specified in the [geom_point()](http://ggplot2.tidyverse.org/reference/geom_point.html) layer (p3).

We can utilize the colour aesthetic to more specifically visualize our data. For example, what if we wanted to know if the 'discovery' or 'extension' cohorts within our data (specified by the 'dataset' variable) had a higher tumor purity? We will use [geom_density()](http://ggplot2.tidyverse.org/reference/geom_density.html) to plot a density kernel of the tumor VAF values, but divide the cohort based upon the dataset subsets (specified within the colour aesthetic).

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

In p1, we told the [geom_density()](http://ggplot2.tidyverse.org/reference/geom_density.html) layer to differentiate the data based upon the 'dataset' column using the colour aesthetic. We see this in our p1 plot, that our result contains two density curves that use two different colored lines to specify our datasets. In p2, we are telling our [geom_density()](http://ggplot2.tidyverse.org/reference/geom_density.html) layer to differentiate the datasets using the "fill" aesthetic instead. We globally assign the line colour ("black") and the fill transparency (alpha=0.75). In addition, we utilize the adjust parameter to reduce the smoothing [geom_density()](http://ggplot2.tidyverse.org/reference/geom_density.html) uses when computing it's estimate. Now, our datasets are specificed by the fill (or filled in color) of each density curve. In p3, we append the [scale_fill_manual()](http://ggplot2.tidyverse.org/reference/scale_manual.html) layer to manually define the fill colours we would like to appear in the plot. This can be done with the colour aesthetic in p1 using [scale_colour_manual()](http://ggplot2.tidyverse.org/reference/scale_manual.html), since we define the colour by the dataset variable.

### Faceting
Depending on the geometric object used there are up to 10 ways to map an aesthetic to a variable. These are with the x-axis, y-axis, fill, colour, shape, alpha, size, labels, and facets. Faceting in ggplot allows us to quickly create multiple related plots at once with a single command. There are two facet commands, [facet_wrap()](http://ggplot2.tidyverse.org/reference/facet_wrap.html) will create a 1 dimensional sequence of panels based on a one sided linear formula. Similarly [facet_grid()](http://ggplot2.tidyverse.org/reference/facet_grid.html) will create a 2 dimensional grid of panels. Let's try and answer a few quick questions about our data using facets.

```R
# what is the most common mutation type among SNP's
p1 <- ggplot(variantData[variantData$type == "SNP",]) + geom_bar(aes(x=trv_type))
p1

# what is the relation of tiers to mutation type
p2 <- ggplot(variantData[variantData$type == "SNP",]) + geom_bar(aes(x=trv_type, fill=tier))
p2

# which reference base is most often mutated
p2 <- p2 + facet_wrap(~reference)

# which transitions and transversions occur most frequently
p2 <- p2 + facet_grid(variant ~ reference)
```

### ggplot themes
Almost every aspect of a ggplot object can be altered. We've already gone over how to alter the display of data but what if we want to alter the display of the non data elements? Fortunately there is a function for that called [theme()](http://ggplot2.tidyverse.org/reference/theme.html). You'll notice in the previous plot some of the x-axis names are colliding with one another, let's fix that and alter some of the theme parameters to make the plot more visually appealing.

```R
# recreate plot p2 if it's not in your environment
p2 <- ggplot(variantData[variantData$type == "SNP",]) + geom_bar(aes(x=trv_type, fill=tier)) + facet_grid(variant ~ reference)
p2

# load in a theme with a few presets set
p3 <- p2 + theme_bw()
p3

# put the x-axis labels on a 45 degree angle
p4 <- p3 + theme(axis.text.x=element_text(angle=45, hjust=1))
p4

# altering a few more visual apects
p5 <- p4 + theme(legend.position="top", strip.text=element_text(colour="white"), strip.background=element_rect(fill="black"))
p5

# Let's remove the y-axis ticks as well
p6 <- p5 + theme(axis.title.x=element_blank())
p6
```

Let's take a few minutes to discuss what is going on here. In p3, we've used [theme_bw()](http://ggplot2.tidyverse.org/reference/ggtheme.html), this function just changes a series of values from the basic default [theme()](http://ggplot2.tidyverse.org/reference/theme.html). There are many such "complete themes" in ggplot and a few external packages as well containing additional "complete themes" such as [ggtheme](https://cran.r-project.org/web/packages/ggthemes/vignettes/ggthemes.html). In p4, we alter the axis.text.x parameter, we can see from the documentation that axis.text.x inherits from [element_text()](http://ggplot2.tidyverse.org/reference/element.html) which is just saying that any parameters in [element_text()](http://ggplot2.tidyverse.org/reference/element.html) also apply to axis.text.x. In this specific case we alter the angle of text to 45 degrees, and set the horizontal justification to the right. In p5 we change the position of the legend, change the colour of the strip.text, and change the strip background. Finally in p6 we remove the x-axis label with [element_blank()](http://ggplot2.tidyverse.org/reference/element.html) which will draw nothing.

### Changing the order of aesthetic mappings
In ggplot the order in which a variable is plotted is determined by the levels of the factor for that variable.
We can view the levels of a column within a dataframe with the [levels()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/levels) command and we can subsequently change the order of those levels with the [factor()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/factor) command. We can then use these to change the order in which aesthetic mappings such as plot facets and discrete axis variables are plotted. Lets look at an example using the previous faceted plots we made (p6)

```R
# recreate plot p6 if it's not in your environment
p6 <- ggplot(variantData[variantData$type == "SNP",]) + geom_bar(aes(x=trv_type, fill=tier)) + facet_grid(variant ~ reference) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", strip.text=element_text(colour="white"), strip.background=element_rect(fill="black"), axis.title.x=element_blank())

# view the order of levels in the reference and trv_type columns
levels(variantData$reference)
levels(variantData$trv_type)

# reverse the order of the levels
variantData$reference <- factor(variantData$reference, levels=rev(levels(variantData$reference)))
variantData$trv_type <- factor(variantData$trv_type, levels=rev(levels(variantData$trv_type)))

# view plot p6 with the new order of variables
p6 <- ggplot(variantData[variantData$type == "SNP",]) + geom_bar(aes(x=trv_type, fill=tier)) + facet_grid(variant ~ reference) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", strip.text=element_text(colour="white"), strip.background=element_rect(fill="black"), axis.title.x=element_blank())
p6
```

We can see that reversing the order of levels in the reference column has subsequently reversed the reference facets (right side of plot). Similarily reversing the order of the trv_type column levels has reversed the labels on the x-axis.

### Saving ggplot2 plots
To save a plot or any graphical object in R, you first have to initalize a graphics device, this can be done with [pdf()](https://www.rdocumentation.org/packages/grDevices/versions/3.4.1/topics/pdf), [png()](https://www.rdocumentation.org/packages/grDevices/versions/3.4.1/topics/png), [svg()](https://www.rdocumentation.org/packages/grDevices/versions/3.4.1/topics/cairo), [tiff()](https://www.rdocumentation.org/packages/grDevices/versions/3.4.1/topics/png), and [jpeg()](https://www.rdocumentation.org/packages/grDevices/versions/3.4.1/topics/png). You can then print the plot as you normally would and close the graphics device using [dev.off()](https://www.rdocumentation.org/packages/grDevices/versions/3.4.1/topics/dev). Alternatively ggplot2 has a function specifically for saving ggplot2 graphical objects called [ggsave()](http://ggplot2.tidyverse.org/reference/ggsave.html). A helpfull tip to keep in mind when saving plots is to allow enough space for the plot to be plotted, if the plot titles for example look truncated try increasing the size of the graphics device. Changes in the aspect ratio between a graphics device height and width will also change the final appearance of plots.

```R
# save the last plot we made.
pdf(file="p1.pdf", height="8", width="11")
p1
dev.off()

# alternatively ggsave will save the last plot made
ggsave("p1.pdf", device="pdf")
```

### ggplot2 Practice examples
Now that we've had an introduction to ggplot2 let's try a few practice example. In the section below we will provide instructions for loading and manipulating a dataset, a plot will then be provided and we ask that you attempt to recreate it. The boxes below will give the answers.

Often it is useful to compare tumor variant allele frequencies among samples to get a sense of the tumor purity and to determine the existense of sub-clonal populations among the tumor. Let's use the [ggplot2ExampleData.tsv](http://www.genomedata.org/gen-viz-workshop/ggplot2ExampleData.tsv) dataset we've been using to explore this.
```R
# load the dataset
variantData <- read.delim("ggplot2ExampleData.tsv")
variantData <- variantData[variantData$dataset == "discovery",]
```
{% include figure.html image="/assets/ggplot2Example1.png" width="950" %}
{% include question.html question="Get a hint!" answer='look at geom_violin(), change labels with xlab() and ylab()'%}
{% include question.html question="What is the code to create the violin plot above?" answer='ggplot() + geom_violin(data=variantData, aes(x=Simple_name, y=tumor_VAF)) + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab("Sample") + ylab("Variant Allele Fraction")'%}

Looking good, but the plot looks dull, try adding some color to the violin plots and let's see where the points for the underlying data actually reside.

{% include figure.html image="/assets/ggplot2Example2.png" width="950" %}
{% include question.html question="Get a hint!" answer='Try using geom_jitter() to offset points'%}
{% include question.html question="What is the code to create the violin plot above?" answer='ggplot(data=variantData, aes(x=Simple_name, y=tumor_VAF)) + geom_violin(aes(fill=Simple_name)) + geom_jitter(width=.1, alpha=.5) + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none") + xlab("Sample") + ylab("Variant Allele Fraction")'%}

Finally let's add some more detail, specifically let's annotate how many points actually make up each violin. The code below will construct the extra data you'll need to make the final plot.

```R
variantDataCount <- count(variantData, "Simple_name")
variantDataMax <- aggregate(data=variantData, tumor_VAF ~ Simple_name, max)
variantDataMerge <- merge(variantDataMax, variantDataCount)
```

{% include figure.html image="/assets/ggplot2Example3.png" width="950" %}
{% include question.html question="Get a hint!" answer='You will need to pass variantDataMerge to geom_text()'%}
{% include question.html question="What is the code to create the violin plot above?" answer='ggplot() + geom_violin(data=variantData, aes(x=Simple_name, y=tumor_VAF, fill=Simple_name)) + geom_jitter(data=variantData, aes(x=Simple_name, y=tumor_VAF), width=.1, alpha=.5) + geom_text(data=variantDataMerge, aes(x=Simple_name, y=tumor_VAF + 5, label=freq)) + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none") + xlab("Sample") + ylab("Variant Allele Fraction")'%}

### Additional resources

* [long vs. wide format](http://seananderson.ca/2013/10/19/reshape.html)
* [ggplot2 wiki](https://github.com/tidyverse/ggplot2/wiki)
* [ggplot2 cheatsheet](https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf)

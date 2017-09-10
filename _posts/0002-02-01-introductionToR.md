---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to R
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-01
---

Over this tutorial series, we will be using R extensively, the language underlying graphical programs [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggvis](https://cran.r-project.org/web/packages/ggvis/index.html), and [GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html). We highly recommend familiarity with R for all students. This page is intended to provide a brief overview of R. However, there are a myriad of resources available that we recommend to supplement the information given here (See Additional Resources below).

## R and Rstudio

{% include figure.html image="/assets/R/Rlogo.svg" position="right" width="250" link="https://www.r-project.org/logo/Rlogo.svg" title="R logo" author="R Foundation" license="CC-BY-SA 4.0" license_link="https://creativecommons.org/licenses/by-sa/4.0/"%}

R is an open source functional programming language developed for statistical computing and graphics. The language includes a number of features that make it ideal for data analysis and visualization and is the primary computing language used in this course.

RStudio in an open source integrated development environment (IDE) written and maintained by RStudio for the R programming language. It is available on all major operating systems and contains a number of features designed to make writing and developing R code more efficient. These features include debugging tools, auto-complete, and syntax highlighting.

## Installation

#### R

To use R go to [https://www.r-project.org/](https://www.r-project.org/), select a mirror via the CRAN link located on the top right, download the appropriate binary distribution for your operating system, and follow the on screen instructions when opening the downloaded file. Once R is installed you can call the function [sessionInfo()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/sessionInfo) to view the version of R and all loaded packages.

#### Rstudio

Rstudio can be downloaded at [https://www.rstudio.com/](https://www.rstudio.com/), select the Products->Rstudio tab within the navigation links at the top, then download the "RStudio Desktop - Open Source Edition", and follow the on screen instructions. Rstudio will automatically detect the current R on your operating system and call the R executable on startup.

## CRAN and Bioconductor

#### CRAN

Part of what makes R an attractive option for data analysis and visualization is the support it receives from a large community of developers. Part of this support comes from users adding additional functionality to core R via functions and data through packages. The Comprehensive R Archive network (CRAN) is a network of servers that stores R, documentation, and many of the packages available for R. To date there are 9,975 [packages on CRAN](https://cran.r-project.org/web/packages/), these packages can be installed by running the install.packages() command in an R terminal. For example, to install the ['plyr'](https://cran.r-project.org/web/packages/plyr/index.html) package, run the following command in an R terminal:
```R
# install the plyr package by Hadley Wickham
install.packages("plyr")
```

#### Bioconductor

Bioconductor is another archive of R packages specific to bioinformatics and genomics. This archive is maintained by the Bioconductor core team and is updated bi-annually. To date, there are 2,541 [packages available via bioconductor](http://bioconductor.org/packages/release/BiocViews.html#___Software). Packages are categorized within biocViews as 'Software', 'AnnotationData', and 'ExperimentData.' You can explore the packages available on this webpage by searching within these tags, categorizing packages based upon relevant topics or types of analysis. Bioconductor packages are managed and installed using the [biocLite()](https://www.rdocumentation.org/packages/BiocInstaller/versions/1.22.3/topics/biocLite) function. Note that this function must be sourced (loaded from source code) before trying to install any Bioconductor packages. Bioconductor packages can be installed by running the following in an R terminal:

```R
# Install core bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite()

# Install specific bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")

# upgrade bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
```
## Documentation and example data in R

As with any software, documentation is key to the usefullness for a user. R makes finding documentation very easy. In order to find the documentation for a specific function, simple enter "?" followed by the function name. This will pull up a manual specific to that function. If you enter this into the Rstudio terminal, the function's documentation immediately appears in the 'Help' tab of the lower right window of the screen. In addition, many packages have additional documentation in the form of vignettes. To view these vignettes from R use the [vignette()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/vignette) function, specifying the package name within the function call (e.g. vignette("grid")). The source code for any function in R can also be viewed by typing the function name into R without any parameters or parentheses.

R also has a variety of datasets pre-loaded. In order to view these, type [data()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/data) in your R terminal. We will be using a few of these data sets to illustrate key concepts in this lesson.

## Assignment and data types

When working in any programming language, values are stored as variables. Defining a variable tells the language to allocate space in memory to store that variable. In R, a variable is assigned with the assignment operator "<-" or "=", which assigns a value to the variable in the user workspace or the current scope respectively. For the purposes of this course, assignment should always occur with the "<-" operator. All variables are stored in objects. The least complex object is the atomic vector. Atomic vectors contain values of only one specific data type. Atomic vectors come in several different flavors (object types) defined by their data type. The six main data types (and corresponding object types) for atomic vectors are: "double (numeric)", "integer", "character", "logical", "raw", and "complex". The data type can be checked with the [typeof()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/typeof) function. The object type can be checked with the [class()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/class) function. Alternatively, data type and object type can also be checked with the is.*() family of functions which will return a logical vector (True/False). Object and data types can also be coerced from one type to another using the as.*() family of functions. An example of each data type is shown below. In each case we will check the object and data type. We will also try coercing some objects/data from one type to another. Understanding and moving between data/object types is important because many R functions expect inputs to be of a certain type and will produce errors or unexpected results if provided with the wrong type.

```R
# numeric
foo <- 1.0
is.numeric(foo)
is.double(foo)
typeof(foo)
class(foo)

# integer values are defined by the "L"
bar <- 1L
is.integer(bar)
typeof(bar)
class(bar)

# character, used to represent strings
baz <- "a"
is.character(baz)
typeof(baz)
class(baz)

# logical values are either TRUE or FALSE
qux <- TRUE
is.logical(qux)
typeof(qux)
class(qux)

# the charToRaw() function is used to store the string in a bit format
corge <- charToRaw("Hello World")
is.raw(corge)
typeof(corge)
class(corge)

# complex
grault <- 4 + 4i
is.complex(grault)
typeof(grault)
class(grault)

#Try coercing foo from a numeric vector with data type double to an integer vector
fooi <- as.integer(foo)
typeof(fooi)
class(fooi)

```

## Data structures (objects)
Data structures in R are objects made up of the data types mentioned above. The type of data structure is dependent upon the homogeneity of the stored data types and the number of dimensions. Commonly used data structures include vectors, lists, matrices, arrays, factors, and dataframes. The most common data structure in R is the vector, which contains data in 1 dimension. There are two types: atomic vectors (discussed above), which contain one data type (i.e. all numeric, all character, etc.) and lists, which contain multiple data types. Atomic vectors are created with the [c()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/c) function. Recall from above that the data type contained within an atomic vector can be determined using the [typeof()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/typeof) function and the type of data object/structure determined with [class()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/class) function. These functions can also be used on more complex data structures. Vectors in R can be sliced (extracting a subset) with brackets [[]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract.data.frame), using either a boolean vector or a numeric index.

```R
# Create an atomic vector
vec <- c(1:10)

# test that it is atomic, of object type numeric, of data type integer, and report the data and object type
is.atomic(vec)
is.numeric(vec)
is.integer(vec)
typeof(vec)
class(vec)

# coerce the numeric vector to character
vec <- as.character(vec)
is.character(vec)

# extract the first element of the vector
vec[1]

# extract the character 5
vec[vec == "5"]

# determine which element of the vector contains a "5"
which(vec == "5")
```

Lists are created using the [list()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/list) function and are used to make more complicated data structures. As mentioned, lists can be heterogeneous, containing multiple data types, objects, or structures (even other lists). Like vectors, items from a list can also be extracted using brackets [[]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract). However, single brackets [[]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract) are used to return an element of the list as a list. Double brackets [[[]]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract) are used to return the the designated element from the list. In general, you should always use double brackets [[[]]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract) when you wish to extract a single item from a list in its expected type. You would use the single brackets [[]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract) if you want to extract a subset of the list.

```R
# create list and verify its data and object type
myList <- list(c(1:10), matrix(1:10, nrow=2), c("foo", "bar"))
typeof(myList)
class(myList)

# extract the first element of the list (the vector of integers)
myList[[1]]
typeof(myList[[1]])
class(myList[[1]])

# extract a subset (e.g., the first and third items) of the list into a new list
myList[c(1,3)]
```

It is important to address attributes in our discussion of data structures. All objects can contain attributes, which hold metadata regarding the object. An example of an attribute for vectors are names. We can give names to each element within a vector with the [names()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/names) function, labeling each element of the list with a corresponding name in addition to its index. Another type of data structure is a factor, which we will use extensively in ggplot2 to define the order of categorical variables. Factors hold metadata (attributes) regarding the order and the expected values. A factor can be specifically defined with the function [factor()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/factor). The expected values and order attributes of a factor are specified with the "levels" param.

```R
# create a named vector
vec <- c("foo"=1, "bar"=2)

# change the names of the vector
names(vec) <- c("first", "second")

# coerce the vector to a factor object
vec <- factor(vec, levels=c(2, 1))
```

## Importing and exporting data

As we have seen, data can be created on the fly in R with the various data structure functions. However it is much more likely that you will need to import data into R to perform analysis. Given the number of packages available, if a common filetype exists (XML, JSON, XLSX) R can probably import it. The most common situation is a simple delimited text file. For this the [read.table()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/read.table) function and its various deriviatives are immensely useful. Type in ?read.table in your terminal for more information about its usage. Once data has been imported and analysis completed, you may need to export data back out of R. Similar to [read.table()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/read.table), R has functions for this purpose. [write.table()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/write.table) will export the data given, in a variety of simple delimited text files, to the file path specified. See ?write.table for more information. Two common issues with importing data use these functions have to do with unrecognized missing/NULL values and unexpected data type coercion by R. Missing values are considered in a special way by R. If you use a notation for missing values that R doesn't expect (e.g., "N/A") in a data column that otherwise contains numbers, R may import that column as a vector of character values instead of a vector of numeric/integer values with missing values. To avoid this, you should always list your missing value notations using the `na.strings` parameter or use the default "NA" that R assumes. Similarly, R will attempt to recognize the data structure type for each column and coerce it to that type. Often with biological data, this results in creation of undesired factors and unexpected results down stream. This can be avoided with the `as.is` parameter.


```R

# import data from a tab-delimited file hosted on the course data server
data <- read.table(file="http://genomedata.org/gen-viz-workshop/ggplot2ExampleData.tsv", header=TRUE, sep="\t", na.strings = c("NA","N/A","na"), as.is=c(1:27,29:30))

# view the first few rows of the imported data
head(data)

# create a new dataframe with just the first 10 rows of data and write that to file
subsetdata <- data[1:10,]
outpath <- getwd()
write.table(x=subsetdata, file=paste(outpath,"/subset.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

```

## Data frames, slicing and manipulation

Within this course, the majority of our analysis will involve analyzing data in the structure of data frames. This is the input ggplot2 expects and is a common and useful data structure throughout the R language. Data frames are 2 dimensional and store vectors, which can be accessed with either single brackets [[]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract.data.frame) or the [$](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract.data.frame) operator. When using single brackets [[]](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract.data.frame), a comma is necessary for specifying rows and columns. This is done by calling the data frame [row(s), column(s)]. The [$](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extract.data.frame) operator is used to specify a column name or variable within the data frame. In the following example, we will use the `mtcars` dataset, one of the preloaded datasets within the R framework. "cyl" is one of the categorical variables within the mtcars data frame, which allows us to specifically call an atomic vector.

Data frames can be created using the [data.frame()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/data.frame) function and is generally the format of data imported into R. We can learn about the format of our data frame with functions such as [colnames()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/row%2Bcolnames), [rownames()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/row%2Bcolnames), [dim()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/dim), [nrow()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/nrow), [ncol()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/nrow), and [str()](https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/str). The example below shows the usefulness of some of these functions, but please use the help documentation for further information. Data frames can be combined in R using the [cbind()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/cbind) and [rbind()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/rbind) functions assuming the data frames being combined have the same column or row names, respectively. If they do not, functions exist within various packages to bind data frames and fill in NA values for columns or rows that do no match (refer to the plyr package for more information).

```R
# view the column names of a dataframe
colnames(mtcars)

# view row names of a dataframe
rownames(mtcars)

# view a summary of the dataframe
str(mtcars)

# subset dataframe to cars with only 8 cyl
mtcars[mtcars$cyl == 8,]

# subset dataframe to first two columns
mtcars[,1:2]
```

## Counting and aggregating data

During the course of an analysis, it is often useful to determine the frequency of an event. A function exists for this purpose in the plyr package called [count()](https://www.rdocumentation.org/packages/plyr/versions/1.8.4/topics/count). Here is an example of how we can apply the [count()](https://www.rdocumentation.org/packages/plyr/versions/1.8.4/topics/count) function to the internal R datasets 'iris' and 'mtcars.'

```R
# first load the plyr package if it's not loaded already
library(plyr)

# How many replicates are there for each species of the iris data?
count(iris$Species)

# How many cars in the mtcars dataset have both 8 cylinders and 4 carburetors?
count(mtcars, c("cyl", "carb"))
```

We can also use the [aggregate()](https://www.rdocumentation.org/packages/stats/versions/3.4.1/topics/aggregate) function to splice our data frames and perform more complicated analyses. The [aggregate()](https://www.rdocumentation.org/packages/stats/versions/3.4.1/topics/aggregate) function requires a formula, by which to splice the data frame, and a function to apply to the subsets described in the formula. In the example below, we will find the average displacement (disp) of cars with 8 cylinders (cyl) and/or 4 carburetors (carb). We will use formulas to splice the data frame by displacement and number of cylinders (disp~cyl) or displacement and number of cylinders and carburetors (disp~cyl + carb) and apply the function [mean()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/mean) to each subset.

```R
# find the mean displacement based on the number of cylinders
aggregate(data=mtcars, disp~cyl, mean)

# find the mean displacement based on the number of cylinders and carburetors
aggregate(data=mtcars, disp~cyl + carb, mean)
```

## Apply family of functions

If you are familiar with other coding languages, you are probably comfortable with looping through the elements of a data structure using functions such as [for()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Control) and [while](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Control). The [apply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/apply) family of functions in R make this much easier (and faster) by inherently looping through a data structure without having to increment throught the indices of the structure. The [apply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/apply) family consists of [lapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply) and [sapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply) for lists, vectors, and data frames and [apply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/apply) for data frames. [lapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply) loops through a list, vector, or data frame, and returns the results of the applied function as a list. [sapply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/lapply) loops through a list, vector, or data frame, and returns the results of the function as a vector or a simplified data structure (in comparison to lapply).

```R
# set a seed for consistency
set.seed(426)

# create a list of distributions
x <- list("dist1"=rnorm(20, sd=5), "dist2"=rnorm(20, sd=10))

# find the standard deviation of each list element
lapply(x, sd)

# return the simplified result
sapply(x, sd)
```

[apply()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/apply) loops over a data frame or matrix, applying the specified function to every row (specified by '1') or column (specified by '2'). For example, given a matrix x, apply(x, 1, min) will apply the [min()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extremes) function to every row of the matrix x, returning the minimum value of each row. Similarly, apply(x, 2, min) will apply the [min()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Extremes) function to every column of matrix x, returning the minimum value in each column.

{% include figure.html image="/assets/R/applyTable.png" width="750" %}

```R
# set a seed for consistency
set.seed(426)

# create a matrix
x <- matrix(runif(n=40, min=1, max=100), ncol=5)

# find the minimum value in each row
apply(x, 1, min)

# find the minimum value in each column
apply(x, 2, min)
```

## Functions in R

A function is a way to store a piece of code so we don't have to type the same code repeatedly. Many existing functions are available in base R or the large number of packages available from CRAN and BioConductor. Occasionally however it may be helpful to define your own custom functions. Combining your own functions with the apply commands above is a powerful way to complete complex or custom analysis on your biological data. For example, suppose we want to determine the number of values in a vector above some cutoff.

```R

# Create a vector of 10 randomly generated values between 1 and 100
vec <- runif(n=10, min=1, max=100)

# Create a function to determine the number of values in a vector greater than some cutoff
myfun <- function(x,cutoff){
 len <- length(which(x>cutoff))
 return(len)
}

# Run your the custom function on the vector with a cutoff of 50
myfun(vec,50)

# Create a matrix of 50 randomly generated values, arranged as 5 columns of 10 values
mat <- matrix(runif(n=50, min=1, max=100), ncol=5)

# Now, use the apply function, together with your custom function
# Determine the number of values above 50 in each row of the matrix
apply(mat, 1, myfun, 50)

```

## Basic control structures

In general, it is preferable to "vectorize" your R code wherever possible. Using functions which take vector arguments will greatly speed up your R code. Under the hood, these "vector" functions are generally written in compiled languages such as C++, FORTRAN, etc. and the corresponding R functions are just wrappers. These functions are much faster when applied over a vector of elements because the compiled code is taking care of the low level processes, such as allocating memory rather than forcing R to do it. There are cases when it is impossible to "vectorize" your code. In such cases there are control structures to help out. Let's take a look at the syntax of a for loop to sum a vector of elements and see how it compares to just running [sum()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/sum).

```R
# install and load a benchmarking package
install.packages("microbenchmark")
library(microbenchmark)

# create a numeric vector to sum
x <- c(2, 4, 6, 8, 10)

# write a function to sum all numbers
mySum <- function(x){

    if(!is.numeric(x)){
        stop("non-numeric argument")
    }

	y <- 0
	for(i in 1:length(x)){
		y <- x[i] + y
	}

	return(y)
}

# both functions produce the correct answer
mySum(x)
sum(x)

# run benchmark tests
microbenchmark(mySum(x), sum(x), times = 1000L)
```

In mySum(), we use a for loop to sum all the elements of the vector. The syntax is fairly straightforward. We loop over the length of the argument passed to x and designate i as the variable to store the iteration of the loop. Prior to that, we use an if statement to make sure the user has supplied only numeric values. This statement simply executes the block of code in curly brackets. If the expression in parenthesis is TRUE, we use an [!](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/Logic) to reverse the outcome of the result given by [is.numeric()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/numeric). All of this is defined as a function. These benchmark tests show that [sum()](https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/sum) is 2-3 orders of magnitude faster than our handwritten mySum() function.

## Additional resources

* [http://www.noamross.net/blog/2014/4/16/vectorization-in-r--why.html](http://www.noamross.net/blog/2014/4/16/vectorization-in-r–why.html)

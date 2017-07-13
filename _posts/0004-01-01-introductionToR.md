---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to R
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-01-01
---

Over this tutorial series we be extensively using R. This is the language which underlies with graphical programs ggplot2, ggvis, and GenVisR. As such a familiarity with R is highly reccommended for students. To this end this topic page intended to provide a brief overview of R however there are a myriad of resources available that we reccomend users investigate to supplement the information given here.

## R and Rstudio

{% include figure.html image="/assets/Rlogo.svg" position="right" width="250" link="https://www.r-project.org/logo/Rlogo.svg" title="R logo" author="R Foundation" license="CC-BY-SA 4.0" license_link="https://creativecommons.org/licenses/by-sa/4.0/"%}

R is an open source functional programming language developed for statistical computing and graphics. The language includes a number of features which make it ideal for data analysis and visualization and is the primary computing language used in this course.

RStudio in an open source integrated development environment (IDE) written and maintained by RStudio for the R programming language. It is available on all major operating systems and contains a number of features designed to make writing and developing R code more efficient. These features include debugging tools, auto-complete, syntax highlighting.

## Installation

#### R

To use R go to [https://www.r-project.org/](https://www.r-project.org/), select a mirror via the CRAN link located on the top right, download the appropriate binary distribution for your operating system, and follow the on screen instructions when opening the downloaded file.

#### Rstudio

Rstudio can be downloaded at [https://www.rstudio.com/](https://www.rstudio.com/), select the Products->Rstudio tab within the navigation links at the top, then download the "RStudio Desktop - Open Source Edition", and follow the on screen instructions. Rstudio will automatically detect the current R on your operating system and call the R executable on startup.

## CRAN and Bioconductor

#### CRAN

Part of what makes R an attractive option for data analysis and visualization is the support it receives from a large community of developers. Part of this support comes from users adding additional functionality to core R via functions and data through packages. The Comprehensive R Archive network (CRAN) is a network of servers that stores R, documentation, and many of the packages available for R. To date there are 9,975 [packages on CRAN](https://cran.r-project.org/web/packages/), these packages can be installed by running the following in an R terminal:
```R
# install the plyr package by Hadley Wickham
install.packages("plyr")
```

#### Bioconductor

Bioconductor is another archive of R packages specific to bioinformatics and genomics. The archive is maintained by the bioconductor core team and is updated bi-annually. To date there are 2,541 [packages available via bioconductor](http://bioconductor.org/packages/release/BiocViews.html#___Software). Packages can be found using biocViews which act as tags for relevant topics. Bioconductor packages can be installed by running the following in an R terminal:

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

As with any software, documentation is key to the usefullness for a user. In this regard R excels, to find the documentation for a specific function in R simply prepend a ? to the function name. This will pull up a manual specific to that function. Further many packages have additional documentation in the form of vignettes, to view these vignettes from R use the vignette() function. In addition the source code for any function in R can be viewed by typing the function name into R without any params or parentheses. In addition to it's exceptional documentation R has a variety of data sets pre-loaded, in order to view these simply type data() in your R terminal. We will be using a few of these data sets to illustrate key concepts in this lesson.

## Assignment and Data types

When working in any language you need to store values into variables. Defining a variable tells the language to allocate space in memory to store that variable. In R A variable is assigned with the assignment operator "<-" or "=", which assign in the user workspace or the current scope respectively (For the purposes of this course assignment should always occur with the "<-" operator). All variables are stored in objects the least complex of which is the atomic vector, these in turn are of a specific type. There are six main data types which are: "numeric", "integer", "character", "logical", "raw", and "complex". The data type can be checked with the is.foo() family of function which will return a logical vector, alternatively the data type can be determined with the class() function. An example of each data type is shown below.

```R
# numeric
foo <- 1.0
is.numeric(foo)
class(foo)

# integer values are defined by the "L"
bar <- 1L
is.integer(bar)

# character, used to represent strings
baz <- "a"
is.character(baz)

# logical values are either TRUE or FALSE
qux <- TRUE
is.logical(qux)

# the charToRaw() function is used to store the string in a bit format
corge <- charToRaw("Hello World")
is.raw(corge)

# complex
grault <- 4 + 4i
is.complex(grault)
```

## Data structures
Data structures in R are objects which hold the data types mentioned above. The type of data structure to use depends on the homogeneity of the data types stored and the number of dimensions needed. The most common data structure in R is the vector which contains data in 1 dimension and comes in two types called Atomic vectors and lists. Atomic vectors are homogeneous in that all the data contained within them must be of the same type (i.e. all numeric, all character, etc.). In contrast lists can contain a mix of data types and can even contain other data structures. Atomic vectors are created with the c() function, the type of atomic vector can be determined with the typeof() function. Vectors in R can be spliced with the [] brackets, using either a boolean vector, or a numeric index.

```R
# Create an atomic vector
vec <- c(1:10)

# test that it is atomic, and of type numeric
is.atomic(vec)
is.numeric(vec)

# coerce the numeric vector to character
vec <- as.character(vec)
is.character(vec)

# extract the first element of the vector
vec[1]

# extract the character 5
vec[vec == "5"]
```

Lists are created using the list() function and are used to make more complicated data structures. As mentioned lists can be heterogeneous in regards to data or object type and can even store other lists. Items from a list can also be extracted using the [] brackets, using a single [] bracket will return an element of the list as a list, using [[]] will return the actual element. In general, when working with list you always want to use the double square brackets.

```R
# create list
myList <- list(c(1:10), matrix(1:10, nrow=2), c("foo", "bar"))

# extract the first element of the list
myList[[1]]
```

It is important to cover attributes in our discussion of data structures. All objects can contain attributes which are used to hold metadata regarding the object. An example of an attribute for vectors are names, we can give names to each element within a vector with the names function. Another attribute is a factor, which is used extensively in ggplot2 to determine order. Factors hold metadata regarding the order and the expected values within a vector. they are defined with the function factor().

```R
# create a named vector
vec <- c("foo"=1, "bar"=2)

# change the names of the vector
names(vec) <- c("first", "second")

# coerce the vector to a factor object
vec <- factor(vec, levels=c(2, 1))
```

## importing and exporting data

As we have seen, we can create data on the fly in R with the various data structure functions. However it is much more likely that you will need to import data into R to perform analysis. Given the number of packages available, if a common filetype exists (XML, JSON, XLSX) R can probably import it. The most common situation is a simple delimited text file. For this the read.table() function and it's various deriviatives are immenseley usefull. Once data has been imported and analysis completed you will need to get data back out of R. Similar to read.table(), R has functions for this purpose. write.table() will export the data given to the file on disk specified.

## Data Frames, slicing and manipulation

Within this course the majority of our work will be with data frames. This is the input ggplot2 expects and is a common and useful object throughout the R language. Dataframes are 2 dimensional and store vectors, these vectors can be accessed with either the square brackets or the $ operator. When using [] brackets a comma is needed to specify whether a column or row extracted (see Data Frames, slicing and manipulation). They are created using the data.frame() function and is usually the format of data when reading into R. Data frames can be combined in R using the cbind() and rbind() functions assuming the data frames being combined have the same columns and rows respectively. If they do not functions exist within various packages to bind data frames and fill in NA values for columns or rows that do no match.

```R
# view the column names of a dataframe
colnames(mtcars)

# view row names of a dataframe
rownames(mtcars)

# subset dataframe to cars with only 8 cyl
mtcars[mtcars$cyl == 8,]

# subset dataframe to first two columns
mtcars[,1:2]
```

## Counting and aggregating

During the course of an analysis it is often usefull to determine the frequency of an event. A useful function exists for this purpose in the plyr package ironically called count(). Let's answer a few questions regarding a few internal R data sets using the count function.

```R
# first load the plyr package if it's not loaded already
library(plyr)

# How many replicates are there for each species of the iris data?
count(iris$Species)

# How many cars in the mtcars dataset have both 8 cylinders and 4 carburetors
count(mtcars, c("cyl", "carb"))
```

As we can see count() is exceptionally usefull but what if we want to do something more complicated like find the average displacement of cars with 8 cylinders and 4 carburetors in the mtcars dataset. Luckily theres a function for that, aggregate() will splice data based on a formula and apply a function across those subsets. Lets go over a few examples and work up to what we want.

```R
# find the mean displacement based on the number of cylinders
aggregate(data=mtcars, disp~cyl, mean)

# find the mean displacement based on the number of cylinders and carburetors
aggregate(data=mtcars, disp~cyl + carb, mean)
```

## Apply family of functions

The apply() functions provide the ablility to loop over differing data structures through various ways. There are there a many of these functions however over this course we will primarily use apply(), lapply(), and mapply(). apply() will apply a function over either the rows or columns of a matrix the determination of which is provided by the second argument to the function call. For example given a matrix x, apply(x, 1, min) will apply the min function to every row of the matrix x. Similarily apply(x, 2, min) will apply the min function to every column of the matrix x.

{% include figure.html image="/assets/applyTable.png" width="750" %}

```R
# set a seed for consistency
set.seed(426)

# create a matrix
x <- matrix(runif(40, 1, 100), ncol=5)

# find the minimum value in each row
apply(x, 1, min)

# find the minimum value in each column
apply(x, 2, min)
```

lapply() loops over either a list or vector of elements and returns the results of the applied funciton as a list. sapply() is similar and is in fact a wrapper for lapply(), the primary difference between the two is that
sapply will simplify the data structure if possible.

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

## Basic control structures

In general it is preferable to "vectorize" your R code wherever possible. Using functions which take vector arguments will greatly speed up your R code. Under the hood these "vector" functions are generally written in compiled languages such as C++, FORTRAN, etc. and the corresponding R functions are just wrappers. These functions are much faster when applied over a vector of elements because the compiled code is taking care of the low level processes such as allocating memory rather than forcing R to do it. There are cases when it is impossible to "vectorize" your code. In such cases there are control structures to help out. Let's take a look at the syntax of a for loop to sum a vector of elements and see how it compares to just running sum().

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

# run benchmark testss
microbenchmark(mySum(x), sum(x))
```

In mySum we use a for loop to sum all the elements of the vector. The syntax is fairly straight forward, we loop over the length of the argument passed to x and designate i as the variable to store the iteration of the loop. Prior to that we use an if statement to make sure the user has supplied a numeric vector. This statement simply executes the block of code in curly brackets if the expression in parenthesis is TRUE, we use an ! to reverse the outcome of the result give by is.numeric(). All of this is defined as a function which is just a way to store a piece of code so we don't have to type the same code over and over. As we can see sum() is 2-3 orders of magnitude faster.

## Additional resources

* [http://www.noamross.net/blog/2014/4/16/vectorization-in-r--why.html](http://www.noamross.net/blog/2014/4/16/vectorization-in-r–why.html)

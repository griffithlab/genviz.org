---
title: Introduction to R
categories:
    - Day 2
feature_image: "https://unsplash.it/1200/400?image=200"
date: 0004-01-01
---

## R and Rstudio

R is an open source functional programming language developed for statistical computing and graphics. The language includes a number of features which make it ideal for data analysis and visualization and is the primary computing language used in this course.

RStudio in an open source integrated development environment (IDE) written and maintained by RStudio for the R programming language. It is available on all major operating systems and contains a number of features designed to make writing and developing R code more efficient. These features include debugging tools, auto-complete, syntax highlighting.

## Installation

### R

To use R go to https://www.r-project.org/, select a mirror via the CRAN link located on the top right, download the appropriate binary distribution for your operating system, and follow the on screen instructions when opening the downloaded file.

### Rstudio

Rstudio can be downloaded at https://www.rstudio.com/, select the Products->Rstudio tab within the navigation links at the top, then download the "RStudio Desktop - Open Source Edition", and follow the on screen instructions. Rstudio will automatically detect the current R on your operating system and call the R executable on startup.

## CRAN and Bioconductor

### CRAN

Part of what makes R an attractive option for data analysis and visualization is the support it receives from a large community of developers. Part of this support comes from users adding additional functionality to core R via functions and data through packages. The Comprehensive R Archive network (CRAN) is a network of servers that stores R, documentation, and many of the packages available for R. To date there are 9,975 [packages on CRAN](https://cran.r-project.org/web/packages/), these packages can be installed by running the following in an R terminal:
```R
# install the plyr package by Hadley Wickham
install.packages("plyr")
```

### Bioconductor

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

## Assignment and Data types

When working in any language you need to store values into variables. Defining a variable tells the language to allocate space in memory to store that variable. In R A variable is assigned with the assignment operator "<-" or "=", which assign in the user workspace or the current scope respectively (For the purposes of this course assignment should always occur with the "<-" operator). All variables are stored in objects the least complex of which is the atomic vector, these in turn are of a specific type. There are six main data types which are: "numeric", "integer", "character", "logical", "raw", and "complex". The data type can be checked with the is.foo() family of function which will return a logical vector, alternatively the data type can be determined with the class() function. An example of each data type is shown below:
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
Data structures in R are objects which hold the data types mentioned above. The type of data structure to use depends on the homogeneity of the data types stored and the number of dimensions needed. The most common data structure in R is the vector which contains data in 1 dimension and comes in two types called Atomic vectors and lists. Atomic vectors are homogeneous in that all the data contained within them must be homogeneous (i.e. all numeric, all character, etc.). In contrast lists can contain a mix of data types and can even contain other data structures. Atomic vectors are created with the c() function, the type of atomic vector can be determined with the typeof() function. Vectors in R can be spliced with the [] brackets, using either a boolean vector, or a numeric index.
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
It is important to cover attributes in our discussion of data structures. All objects can contain attributes which are used to hold metadata regarding the object. An example of an attribute for vectors are names, we can give names to each element within a vector with the names function. Another attribute is a factor, which is used extensively in ggplot2 to determine order. Factors hold metadata regarding the order and the expected values within a vector. they are defined with the factor function()
```R
# create a named vector
vec <- c("foo"=1, "bar"=2)

# change the names of the vector
names(vec) <- c("first", "second")

# coerce the vector to a factor object
vec <- factor(vec, levels=c(2, 1))
```
Within this course the majority of our work will be with data frames. This is the input ggplot2 expects and is common and useful object throughout the R language. Dataframes are 2 dimensional and store vectors, these vectors can be accessed with either the square brackets or the $ function. When using [] brackets a comma is needed to specify whether a column or row extracted (see Data Frames, slicing and manipulation). They are created using the data.frame() function and is usually the format of data when reading into R. Data frames can be combined in R using the cbind() and rbind() functions assuming the data frames being combined have the same columns and rows respectively. If they do not functions exist within the plyr package cbind.fill() and rbind.fill() to bind data frames and fill in NA values for columns or rows that do no match.
## Reading and writing data
## Data Frames, slicing and manipulation
## Counting, aggregating
## Basic control structures
## Apply family of functions
## Additional resources

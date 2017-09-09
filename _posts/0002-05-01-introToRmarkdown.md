---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to Rmarkdown
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-01
---

A usefull feature within the R ecosystem are RMarkdown fies. These files allow a user to intersperse notes with code providing a useful framework for sharing scientific reports in a transparent and easily reproduceable way. They combine both markdown syntax for styling notes and Rscripts for running code and producing figures. Reports can be output in a variety of file formats including HTML documents, PDF documents, and even presentation slides.

### Installing Rmarkdown
To start let's make a simple Rmarkdown file with Rstudio,  you will need to install the [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html) package from cran.
```R
# install rmarkdown
install.packages("rmarkdown")
```

Once [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html) has been installed we can select File -> New File -> R Markdown to set up an intial rmarkdown template.

{% include figure.html image="/assets/Rmarkdown/intro_2_rmarkdown_1.png" width="550" %}

Rstudio will ask you to choose an output format, and to add a title and author, for now we will just use the default HTML format however this can be changed at any time within the Rmarkdown template. Go ahead and select okay when you have added your name and a title.

{% include figure.html image="/assets/Rmarkdown/intro_2_rmarkdown_2.png" width="450" %}

Rsutdio should now have made a template for us, let's go over a few introductory topics in regards to this template. At the top of the file you will see what looks like a YAML header denoted by `---`. This is where the defaults for building the file are set. You will notice that [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html) has pre-populated some fields based on what we supplied when we initalized the markdown template. You can output the rmarkdown document using the `Knit` button in the top left hand corner. This is the same as calling the function `render()` which takes the path to the rmarkdown file as input. This file should end in a .Rmd extension to denote it as an rmarkdown file, though Rstudio will take care of this for you the first time you hit `Knit`. Rstudio also has a convenient way to insert code using the `insert` button to the right. You might notice that not only does [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html) support R, but also bash, python and a few other languages as well. Though in order to work these languages will need to be installed before using `Knit`.

{% include figure.html image="/assets/Rmarkdown/intro_2_rmarkdown_3.png" width="650" %}

As we have mentioned you can insert a code chunk using the `insert` button on the top right. For example as shown below when selecting insert -> R, we get a code chunk formatted for R. However you can also add parameters to this code chunk to alter it's default behavior. A full list of these parameters is available [here](https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf).

{% include figure.html image="/assets/Rmarkdown/intro_2_rmarkdown_4.png" width="650" %}

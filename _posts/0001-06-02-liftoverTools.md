---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to liftover tools
categories:
    - Module-01-Intro
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-06-02
---

A common analysis task is to convert genomic coordinates between different assemblies. Probably the most common situation is that you have some coordinates for a particular version of a reference genome and you want to determine the corresponding coordinates on a different version of the reference genome for that species. For example, you have a bed file with exon coordinates for human build GRC37 (hg19) and wish to update to GRCh38. Many resources exist for performing this and other related tasks. In this section we will go over a few tools to perform this type of analysis, in many cases these tools can be used interchangeably. This post is inspired by this [BioStars post](https://www.biostars.org/p/65558/) (also created by the authors of this workshop).

### Reference assemblies
First let's go over what a reference assembly actually is, in essence it's just a representation of the nucleotide sequence from a cohort. These assemblies allow for a shortcut when mapping reads as they can be mapped to the assembly, rather than each other, to piece the genome of an individual together. This has a number of benefits, the most obvious of which is that it is far more effecient than attempting to build a genome from scratch. By its very nature however using this approach means there is no perfect reference assembly for an individual due to polymorphism (i.e. snps, hla-type, etc.). Further due to the presence of repetitive structural elements such as duplications, inverted repeats, tandem repeats, etc. a given assembly is almost always incomplete, and is constantly being improved upon. This leads to the publication of new assembly versions every so often such as grch37 (Feb. 2009) and grch38 (Dec. 2013). It is also good to be aware that different organization can publish different reference assemblies, for example grch37 (NCBI) and hg19 (UCSC) are identical save for a few minor differences such as in the mitochondria sequence and annotation of chromosomes (1 vs chr1). For a nice summary of genome versions and their release names refer to the [Assembly Releases and Versions FAQ](http://genome.ucsc.edu/FAQ/FAQreleases.html).

### Resources available
There are many resources available to convert coordinates from one assemlby to another, we will go over a few of these however below you will find a more complete list. The [UCSC liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) tool is one of the more popular however choosing one of these will mostly come down to personal preference.

* **UCSC liftOver:** This tool is available through a simple [web interface](http://genome.ucsc.edu/cgi-bin/hgLiftOver) or it can be downloaded as a [standalone executable](http://hgdownload.cse.ucsc.edu/admin/exe/). To use the executable you will also need to download the appropriate [chain file](http://hgdownload.cse.ucsc.edu/downloads.html#liftover). Each chain file describes conversions between a pair of genome assemblies. Liftover can be used through [Galaxy](https://usegalaxy.org/) as well. There is a python implementation of liftover called [pyliftover](https://pypi.python.org/pypi/pyliftover) that does conversion of point coordinates only.

* **NCBI Remap:** This tool is conceptually similar to liftOver in that in manages conversions between a pair of genome assemblies but it uses different methods to achieve these mappings. It is also available through a simple [web interface](http://www.ncbi.nlm.nih.gov/genome/tools/remap) or you can use the [API for NCBI Remap](http://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/api).

* **The Ensembl API:** The final example I described above (converting between coordinate systems within a single genome assembly) can be accomplished with the [Ensembl core API](http://ensembl.org/info/docs/api/core/index.html#api). Many examples are provided within the [installation](http://ensembl.org/info/docs/api/api_installation.html), [overview](http://ensembl.org/info/docs/api/core/core_API_diagram.html), [tutorial](http://ensembl.org/info/docs/api/core/core_tutorial.html) and [documentation](http://ensembl.org/info/docs/Doxygen/core-api/index.html) sections of the Ensembl API project. In particular, refer to these sections of the [tutorial](http://ensembl.org/info/docs/api/core/core_tutorial.html): 'Coordinates', 'Coordinate systems', 'Transform', and 'Transfer'.

* **Assembly Converter:** Ensembl also offers their own simple web interface for coordinate conversions called the [Assembly Converter](https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core).

* **rtracklayer:** For [R](http://cran.us.r-project.org/index.html) users, [Bioconductor](http://www.bioconductor.org/) has an implementation of UCSC liftOver in the [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package. To see documentation on how to use it, open an R session and run the following commands.

* **CrossMap:** A standalone [open source program](http://crossmap.sourceforge.net/) for convenient conversion of genome coordinates (or annotation files) between different assemblies. It supports most commonly used file formats including SAM/BAM, Wiggle/BigWig, BED, GFF/GTF, VCF. CrossMap is designed to liftover genome coordinates between assemblies. It’s not a program for aligning sequences to reference genome. Not recommended for converting genome coordinates between species.

* **Flo:** A [liftover pipeline](https://github.com/yeban/flo) for different reference genome builds of the same species. It describes the process as follows: "align the new assembly with the old one, process the alignment data to define how a coordinate or coordinate range on the old assembly should be transformed to the new assembly, transform the coordinates."

### Using UCSC liftOver
[Sex linkage](https://en.wikipedia.org/wiki/Sex_linkage) was first discovered by [Thomas Hunt Morgan](https://en.wikipedia.org/wiki/Thomas_Hunt_Morgan) in 1910 when he observed that the eye color of **Drosophila melanogaster** did not follow typical mendelian inheritance. This was discovered to be caused by the **white** gene located on chromosome X at coordinates 2684762-2687041 for assembly dm3. Let's use [UCSC liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to determine where this gene is located on the latest reference assembly for this species, dm6. First navigate to the liftOver site at  [https://genome.ucsc.edu/cgi-bin/hgLiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) and set both the original and new genomes to the appropriate species, "D. melanogaster".

{% include figure.html image="/assets/liftOver/liftOver_1.png" width=650 %}

We want to transfer our coordinates from the dm3 assembly to the dm6 assembly so let's make sure the original and new assemblies are set appropriately as well.

{% include figure.html image="/assets/liftOver/liftOver_2.png" width=650 %}

Finally we can paste our coordinates to transfer or upload them in bed format.

{% include figure.html image="/assets/liftOver/liftOver_3.png" width=650 %}

The page will refresh and a results section will appear where we can download the transferred cordinates in bed format.

{% include figure.html image="/assets/liftOver/liftOver_4.png" width=650 %}

Try and compare the old and new coordinates in the [UCSC genome browser]() for their respective assemblies, do they match the same gene?

{% include question.html question="Answer" answer='Yes, both coordinates match the coding sequence for the <i>w</i> gene from transcript CG2759-RA'%}

### Using rtracklayer
In another situation you may have coordinates of a gene and wish to determine the corresponding coordinates in another species. This is a common situation in evolutionary biology where you will need to find coordinates for a conserved gene across species to perform a phylogenetic analysis. Let's use the rtracklayer package on bioconductor to find the coordinates of the **H3F3A** gene located at chr1:226061851-226071523 on the hg38 human assembly in the canFam3 assembly of the canine genome. To start install the [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package from bioconductor, as mentioned this is an R implementation of the UCSC liftover.

```R
# install and load rtracklayer
# source("https://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
library("rtracklayer")
```

The function we will be using from this package is [liftover()](https://www.rdocumentation.org/packages/rtracklayer/versions/1.32.1/topics/liftOver) and takes two arguments as input. The first of these is a [GRanges](https://www.rdocumentation.org/packages/GenomicRanges/versions/1.24.1/topics/GRanges-class) object specifying coordinates to perform the query on. This class is from the [GenomicRanges](https://www.rdocumentation.org/packages/GenomicRanges/versions/1.24.1) package maintained by bioconductor and was loaded automatically when we loaded the [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html) library. The second item we need is a [chain file](https://genome.ucsc.edu/goldenpath/help/chain.html), which is a format which describes pairwise alignments between sequences allowing for gaps. The UCSC website maintains a selection of these on it's [genome data page](http://hgdownload.soe.ucsc.edu/downloads.html). Navigate to this page and select "liftOver files" under the hg38 human genome, then download and extract the "hg38ToCanFam3.over.chain.gz" chain file.

{% include figure.html image="/assets/liftOver/liftOver_5.png" width=650 %}

Next all we need to do is to create our [GRanges](https://www.rdocumentation.org/packages/GenomicRanges/versions/1.24.1/topics/GRanges-class) object to contain the coordinates chr1:226061851-226071523 and import our chain file with the function [import.chain()]. We can then supply these two parameters to [liftover()](https://www.rdocumentation.org/packages/rtracklayer/versions/1.32.1/topics/liftOver).

```R
# specify coordinates to liftover
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=226061851, end=226071523))

# import the chain file
chainObject <- import.chain("hg38ToCanFam3.over.chain")

# run liftOver
results <- as.data.frame(liftOver(grObject, chainObject))
```

How many different regions in the canine genome match the human region we specified?

{% include question.html question="Answer" answer='210, these return the ranges mapped for the corresponding input element' %}

### Exercises
Try to perform the same task we just complete with the web version of [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver), how are the results different?
{% include question.html question="Answer" answer='both methods provide the same overall range, however using rtracklayer is not simplified and contains multiple ranges corresponding to the chain file.' %}

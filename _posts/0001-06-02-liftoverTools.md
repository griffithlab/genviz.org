---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to liftover tools
categories:
    - Module 1
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-06-02
---

A common analysis task is to convert is to convert genomic coordinates between different assemblies. Probably the most common situation is that you have some coordinates for a particular version of a reference genome and you want to determine the corresponding coordinates on a different version of the reference genome for that species. For example, you have a bed file with exon coordinates for human build GRC37 (hg19) and wish to update to GRCh38. Many resources exist for performing this and other related tasks. In this section we will go over a few tools to perform this type of analysis, in many cases these tools can be used interchangeably.

### Reference assemblies
First let's go over what a reference assembly actually is, in essence it's just a representation of the nucleotide sequence from a cohort. These assemblies allow for a shortcut when mapping reads as they can be mapped to the assembly, rather than each other, to piece the genome of an individual together. This has a number of benefits, the most obvious of which is that it is far more effecient than attempting to build a genome from scratch. By its very nature however using this approach means there is no perfect reference assembly for an individual due to polymorphism (i.e. snps, hla-type, etc.). Further due to the presence of repetitive structural elements such as duplications, inverted repeats, tandem repeats, etc. a given assembly is almost always incomplete, and is constantly being improved upon. This leads to the publication of new assembly versions every so often such as grch37 (Feb. 2009) and grch38 (Dec. 2013). It is also good to be aware that different organization can publish different reference assemblies, for example grch37 (NCBI) and hg19 (UCSC) are identical save for a few minor differences such as in the mitochondria sequence and annotation of chromosomes (1 vs chr1).

### Resources available
There are many resources available to accomplish these tasks, we will go over a few of these however below you will find a more complete list. The [UCSC liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) tool is one of the more popular however choosing one will mostly come down to personal preference.

* **UCSC liftOver:** This tool is available through a simple [web interface](http://genome.ucsc.edu/cgi-bin/hgLiftOver) or it can be downloaded as a [standalone executable](http://hgdownload.cse.ucsc.edu/admin/exe/). To use the executable you will also need to download the appropriate [chain file](http://hgdownload.cse.ucsc.edu/downloads.html#liftover). Each chain file describes conversions between a pair of genome assemblies. Liftover can be used through [Galaxy](https://usegalaxy.org/) as well. There is a python implementation of liftover called [pyliftover](https://pypi.python.org/pypi/pyliftover) that does conversion of point coordinates only.

* **NCBI Remap:** This tool is conceptually similar to liftOver in that in manages conversions between a pair of genome assemblies but it uses different methods to achieve these mappings. It is also available through a simple [web interface](http://www.ncbi.nlm.nih.gov/genome/tools/remap) or you can use the [API for NCBI Remap](http://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/api).

* **The Ensembl API:** The final example I described above (converting between coordinate systems within a single genome assembly) can be accomplished with the [Ensembl core API](http://ensembl.org/info/docs/api/core/index.html#api). Many examples are provided within the [installation](http://ensembl.org/info/docs/api/api_installation.html), [overview](http://ensembl.org/info/docs/api/core/core_API_diagram.html), [tutorial](http://ensembl.org/info/docs/api/core/core_tutorial.html) and [documentation](http://ensembl.org/info/docs/Doxygen/core-api/index.html) sections of the Ensembl API project. In particular, refer to these sections of the [tutorial](http://ensembl.org/info/docs/api/core/core_tutorial.html): 'Coordinates', 'Coordinate systems', 'Transform', and 'Transfer'.

* **Assembly Converter:** Ensembl also offers their own simple web interface for coordinate conversions called the [Assembly Converter](https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core).

* **rtracklayer:** For [R](http://cran.us.r-project.org/index.html) users, [Bioconductor](http://www.bioconductor.org/) has an implementation of UCSC liftOver in the [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package. To see documentation on how to use it, open an R session and run the following commands.

* **CrossMap:** A standalone [open source program](http://crossmap.sourceforge.net/) for convenient conversion of genome coordinates (or annotation files) between different assemblies. It supports most commonly used file formats including SAM/BAM, Wiggle/BigWig, BED, GFF/GTF, VCF. CrossMap is designed to liftover genome coordinates between assemblies. It’s not a program for aligning sequences to reference genome. Not recommended for converting genome coordinates between species.

* **Flo:** A [liftover pipeline](https://github.com/yeban/flo) for different reference genome builds of the same species. It describes the process as follows: "align the new assembly with the old one, process the alignment data to define how a coordinate or coordinate range on the old assembly should be transformed to the new assembly, transform the coordinates."

This comes up in several contexts. Probably the most common is that you have some coordinates for a particular version of a reference genome and you want to determine the corresponding coordinates on a different version of the reference genome for that species. For example, you have a bed file with exon coordinates for human build GRC37 (hg19) and wish to update to GRCh38. By the way, for a nice summary of genome versions and their release names refer to the Assembly Releases and Versions FAQ

Or perhaps you have coordinates of a gene and wish to determine the corresponding coordinates in another species. For example, you have coordinates of a gene in human GRCh38 and wish to determine corresponding coordinates in mouse mm10.

Finally you may wish to convert coordinates between coordinate systems within a single assembly. For example, you have the coordinates of a series of exons and you want to determine the position of these exons with respect to the transcript, gene, contig, or entire chromosome.

There are now several well known tools that can help you with these kinds of tasks:

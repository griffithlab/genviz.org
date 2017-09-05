---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to liftover tools
categories:
    - Module 1
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-06-02
---


This comes up in several contexts. Probably the most common is that you have some coordinates for a particular version of a reference genome and you want to determine the corresponding coordinates on a different version of the reference genome for that species. For example, you have a bed file with exon coordinates for human build GRC37 (hg19) and wish to update to GRCh38. By the way, for a nice summary of genome versions and their release names refer to the Assembly Releases and Versions FAQ

Or perhaps you have coordinates of a gene and wish to determine the corresponding coordinates in another species. For example, you have coordinates of a gene in human GRCh38 and wish to determine corresponding coordinates in mouse mm10.

Finally you may wish to convert coordinates between coordinate systems within a single assembly. For example, you have the coordinates of a series of exons and you want to determine the position of these exons with respect to the transcript, gene, contig, or entire chromosome.

There are now several well known tools that can help you with these kinds of tasks:

UCSC liftOver. This tool is available through a simple web interface or it can be downloaded as a standalone executable. To use the executable you will also need to download the appropriate chain file. Each chain file describes conversions between a pair of genome assemblies. Liftover can be used through Galaxy as well. There is a python implementation of liftover called pyliftover that does conversion of point coordinates only.

NCBI Remap. This tool is conceptually similar to liftOver in that in manages conversions between a pair of genome assemblies but it uses different methods to achieve these mappings. It is also available through a simple web interface or you can use the API for NCBI Remap.

The Ensembl API. The final example I described above (converting between coordinate systems within a single genome assembly) can be accomplished with the Ensembl core API. Many examples are provided within the installation, overview, tutorial and documentation sections of the Ensembl API project. In particular, refer to these sections of the tutorial: 'Coordinates', 'Coordinate systems', 'Transform', and 'Transfer'.

Assembly Converter. Ensembl also offers their own simple web interface for coordinate conversions called the Assembly Converter.

Bioconductor rtracklayer package. For R users, Bioconductor has an implementation of UCSC liftOver in the rtracklayer package. To see documentation on how to use it, open an R session and run the following commands.

source("http://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library(rtracklayer)
?liftOver

CrossMap. A standalone open source program for convenient conversion of genome coordinates (or annotation files) between different assemblies. It supports most commonly used file formats including SAM/BAM, Wiggle/BigWig, BED, GFF/GTF, VCF. CrossMap is designed to liftover genome coordinates between assemblies. It’s not a program for aligning sequences to reference genome. Not recommended for converting genome coordinates between species.

Flo. A liftover pipeline for different reference genome builds of the same species. It describes the process as follows: "align the new assembly with the old one, process the alignment data to define how a coordinate or coordinate range on the old assembly should be transformed to the new assembly, transform the coordinates."

---
feature_text: |
  ## Genomic Visualization and Interpretations
title: ggplot2 Continued*
categories:
    - Module-02-R
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-02
---

In the previous section we covered the core aspects of ggplot2, In this section we will cover some of the lesser used features that arre still quite usefull. In addition this section contains additional exercises to reinforce concepts.

To start things off let's go ahead and load in an transcripts annotation database to work with. Bioconductor maintains many of these databases for different species/assemblies, here we load in one from Bioconductor for Hsapiens/HG38. You can view the many different transcript annotation databases bioconductor offers by looking for the [TxDb BiocView](https://bioconductor.org/packages/release/BiocViews.html#___TxDb) on bioconductor.

```R
# install if not already
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# load the library
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
```

It is beyond the purpose of this workshop to explain everything a TxDb object stores. In the code below we are just extract gene related data from chromosome 1-22,X,Y to plot. Note that in the dataframe we produce below the Gene ID is an Entrez ID.

```R
# obtain gene locations from the TxDb object as a data frame
genes <- genes(TxDb, columns=c("gene_id"))
genes <- as.data.frame(genes)
chromosomes <- paste0("chr", c(1:22,"X","Y"))
genes <- genes[genes$seqnames %in% chromosomes,]
```

Okay, we have the core data we'll be using for this section, let's go ahead and use what we know to answer a couple biologically relevant questions. To start let's look at the density of genes across all chromosomes. Try to recreate the plot below.

```R
ggplot(data=genes, aes(x=start + .5*width,)) + geom_histogram(bins=100) + geom_rug(data=genes[genes$strand == "-",], aes(x=start + .5*width), color="tomato3", sides="t", alpha=.1, length=unit(0.1, "npc")) + geom_rug(data=genes[genes$strand == "+",], aes(x=start + .5*width), color="darkorchid4", sides="b", alpha=.1, length=unit(0.1, "npc")) + facet_wrap(~seqnames, scales="free_x") + theme_bw()# + theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
```

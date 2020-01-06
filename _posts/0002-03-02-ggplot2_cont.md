---
feature_text: |
  ## Genomic Visualization and Interpretations
title: ggplot2 Continued*
categories:
    - Module-02-R
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-02
---

In the previous section we covered the core aspects of ggplot2, In this section we provide some additional exercises to reinforce concepts.

To start things off let's go ahead and load in a transcripts annotation database to work with. Bioconductor maintains many of these databases for different species/assemblies, here we load in one from Bioconductor for Hsapiens/HG38. You can view the many different transcript annotation databases bioconductor offers by looking for the [TxDb BiocView](https://bioconductor.org/packages/release/BiocViews.html#___TxDb) on bioconductor.

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

Okay, we have the core data we'll be using for this section, let's go ahead and use what we know to answer a couple biologically relevant questions adding layers as we go to create a more compicated plot. To start let's look at the density of genes across genomic coordinates. Try to recreate the plot below by plotting the density of the genes acorss the genome. Use the gene center coordinate of the gene for this.

{% include question.html question="Get a hint!" answer='To plot the center coordinate of the gene you will need to make another columns in the data frame, you can do this within ggplot or just add another column.'%}
{% include question.html question="What is the code to produce the plot below?" answer='ggplot(data=genes, aes(x=start + .5*width,)) + geom_density()'%}

{% include figure.html image="/assets/ggplot2/ggplot2_cont_density_part1.png" width="950" %}

Pretty easy right? As you would expect gene density is higher twoard the beginning of chromosomes simply because there is more overlap of genomic coordiantes between chromosomes at the star (i.e. all chromosomes start at 1, but are of different lengths). Now let's add to our plot by creating adding a rug layer for genes on the anti-sense strand. Place this rug on the top of the plot, alter the transparency, color, and size of the rug. when your done your plot should look similar to the one below.

{% include question.html question="Get a Hint!" answer='Look at the ggplot2 documentation for geom_rug()'%}
{% include question.html question="What is the code to produce the plot below" answer='ggplot(data=genes, aes(x=start + .5*width,)) + geom_density() + geom_rug(data=genes[genes$strand == "-",], aes(x=start + .5*width), color="tomato3", sides="t", alpha=.1, length=unit(0.1, "npc"))
'%}
{% include figure.html image="/assets/ggplot2/ggplot2_cont_density_part2.png" width="950" %}

While were at it, let's go ahead and add a layer for the sense strand as well, this will go on the bottom of the plot instead of the plot. Go ahead and try and mimic the plot below.

{% include question.html question="What is the code to produce the plot below" answer='ggplot(data=genes, aes(x=start + .5*width,)) + geom_density() + geom_rug(data=genes[genes$strand == "-",], aes(x=start + .5*width), color="tomato3", sides="t", alpha=.1, length=unit(0.1, "npc")) + geom_rug(data=genes[genes$strand == "+",], aes(x=start + .5*width), color="darkorchid4", sides="b", alpha=.1, length=unit(0.1, "npc"))
'%}
{% include figure.html image="/assets/ggplot2/ggplot2_cont_density_part3.png" width="950" %}

We have are basic plot now, but it's still not very informative as all the data from the different chromosomes are colliding with each other. Let's go ahead and fix this by making multiple plots from the same data split out by chromosome. Recall from the privous section that there is a very easy way to do this. Also pay particular attention to how the axis are set for each individual plot and produce the result below.

{% include question.html question="Get a Hint!" answer='look at the ggplot2 documentation for facet_wrap()'%}
{% include question.html question="What is the code to produce the plot below" answer='ggplot(data=genes, aes(x=start + .5*width,)) + geom_density() +
  geom_rug(data=genes[genes$strand == "-",], aes(x=start + .5*width), color="tomato3", sides="t", alpha=.1, length=unit(0.1, "npc")) +
  geom_rug(data=genes[genes$strand == "+",], aes(x=start + .5*width), color="darkorchid4", sides="b", alpha=.1, length=unit(0.1, "npc")) +
  facet_wrap(~seqnames, scales="free")
'%}
{% include figure.html image="/assets/ggplot2/ggplot2_cont_density_part4.png" width="950" %}

We can start to see some interesting trends now, specifically chr6 appeears to have many genes on the p-arm of the chromosome compared to the q-arm. We can also see a couple regions where strand bias might be present, such as in the beginning of chromosome 15. Let's go ahead and finish things up by altering some of the theme aspects of the plot. Reproduce the plot below using theme()

{% include question.html question="What is the code to produce the plot below" answer='ggplot(data=genes, aes(x=start + .5*width,)) + geom_density() + geom_rug(data=genes[genes$strand == "-",], aes(x=start + .5*width), color="tomato3", sides="t", alpha=.1, length=unit(0.1, "npc")) + geom_rug(data=genes[genes$strand == "+",], aes(x=start + .5*width), color="darkorchid4", sides="b", alpha=.1, length=unit(0.1, "npc")) + facet_wrap(~seqnames, scales="free") + theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
'%}
{% include figure.html image="/assets/ggplot2/ggplot2_cont_density_part5.png" width="950" %}

Let's try another exercise, what if we don't care about the density of genes across chromosomes, but instead just want to know which chromosome has the highest gene burden. The code below will produce the data ggplot2 will need to plot this information

```R
geneFreq <- plyr::count(genes$seqnames)
names(geneFreq) <- c("seqnames", "freq")
chrLength <- as.data.frame(seqlengths(TxDb))
chrLength$seqnames <- rownames(chrLength)
colnames(chrLength) <- c("length","seqnames")
chrGeneBurden <- merge(geneFreq, chrLength, all.x=TRUE, by="seqnames")
chrGeneBurden$gene_per_mb <- chrGeneBurden$freq/chrGeneBurden$length * 1000000
```

Go ahead and try and reproduce the plot below, remember to keep adding layers to get closer to the final product and refer to the ggplot2 documention for help.

{% include question.html question="What is the code to produce the plot below" answer='myChrOrder <- as.character(chrGeneBurden[order(chrGeneBurden$gene_per_mb),]$seqnames)
chrGeneBurden$seqnames <- factor(chrGeneBurden$seqnames, levels=myChrOrder)<br><br>ggplot(chrGeneBurden, aes(seqnames, gene_per_mb)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = median(chrGeneBurden$gene_per_mb), linetype=2, color="tomato1") +
  geom_text(aes(label=freq), angle=-55, hjust=1, nudge_y=.5) +
  scale_x_discrete(expand=c(.05, .05)) +
  scale_y_continuous(sec.axis = dup_axis(name="")) +
  geom_point(aes(y=30, x=1), alpha=0) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ylab("Genes Per MB")
'%}
{% include figure.html image="/assets/ggplot2/ggplot2_cont_density_part6.png" width="950" %}

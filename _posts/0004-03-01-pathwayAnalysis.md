---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Pathway analysis
categories:
    - Module-04-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-03-01
---

In the previous section we examined differential expression of genes from the [E-GEOD-50760](https://www.ncbi.nlm.nih.gov/pubmed/25049118) data set. In this section we will use the [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) package to determine if there are any coordinated differential expression patterns in the data set we used for differential expression,  [E-GEOD-50760](https://www.ncbi.nlm.nih.gov/pubmed/25049118).

### What is gage?
generally applicable gene-set enrichment ([gage](https://bioconductor.org/packages/release/bioc/html/gage.html)) is a popular bioconductor package for performing gene-set and pathway analysis. The package works independent of sample sizes, experimental designs, assay platforms, and is applicable to both microarray and rnaseq data sets. In this section we will use [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) and gene sets from the "Kyoto Encyclopedia of Genes and Genomes" ([KEGG](http://www.kegg.jp/)) and "Gene Ontology" ([GO](http://www.geneontology.org/)) databases to perform pathway analysis. Let's go ahead and install [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) and load the differential expression results from the previous section.

```R
# install gage
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("gage", version = "3.8")
library(gage)

# load the differential expression results fro the previous section
load(url("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/deseq2Data_v1.RData"))

# extract the results from the deseq2 data
library(DESeq2)
tumor_v_normal_DE <- results(deseq2Data, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))
```

### setting up gene set databases
In order to perform our pathway analysis we need a list of pathways and their respective genes. The most common databases for this type of data are [KEGG](http://www.kegg.jp/) and [GO](http://www.geneontology.org/). The [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) package has two functions for querying this information in real time, [kegg.gsets()](https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/kegg.gsets) and [go.gsets()](https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/go.gsets), both of which take a species as an argument and will return a list of gene sets and some helpful meta information for subsetting these list. For the [KEGG](http://www.kegg.jp/) database object `kg.hsa$kg.sets` stores all gene sets for the queried species; `kg.hsa$sigmet.idx` and `kg.hsa$dise.idx` store the indices for those gene sets which are classified as signaling and metabolism and disease respectively. We use this information to extract a list of gene sets for the signaling and metabolism and disease subsets. A similar process is used for the [GO](http://www.geneontology.org/) gene sets splitting the master gene set into the three gene ontologies: "Biological Process", "Molecular Function", and "Cellular Component".

```R
# set up kegg database
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]
```

### annotating genes
We have our gene sets now however if you look at one of these objects containing the gene sets you'll notice that each gene set contains a series of integers. These integers are actually [entrez](https://www.ncbi.nlm.nih.gov/gquery/) gene identifiers which presents a problem as our DESeq2 results use ensemble ID's as gene identifiers. We will need to convert our gene identifiers to the same format before we perform the pathway analysis. Fortunately bioconductor maintains genome wide annotation data for many species, you can view these species with the [OrgDb bioc view](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb). This makes converting the gene identifiers relatively straight forward, below we use the [mapIds()](https://www.rdocumentation.org/packages/OrganismDbi/versions/1.14.1/topics/MultiDb-class) function to query the OrganismDb object for the gene symbol, entrez id, and gene name based on the ensembl id. Because there might not be a one to one relationship with these identifiers we also use `multiVals="first"` to specify that only the first identifier should be returned in such cases.

```R
# load in libraries to annotate data
source("https://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi","org.Hs.eg.db"))
library(AnnotationDbi)
library(org.Hs.eg.db)

# annotate the deseq2 results with additional gene identifiers
tumor_v_normal_DE$symbol <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
tumor_v_normal_DE$entrez <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
tumor_v_normal_DE$name <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="GENENAME", keytype="ENSEMBL", multiVals="first")
```

### Preparing DESeq2 results for gage
Before we perform the actuall pathway analysis we need to format our differential expression results into a format suitable for the [gage]() package. Basically this means obtaining the normalized log2 expression values and assigning entrez gene identifiers to these values.

```R
# grab the log fold changes for everything
tumor_v_normal_DE.fc <- tumor_v_normal_DE$log2FoldChange
names(tumor_v_normal_DE.fc) <- tumor_v_normal_DE$entrez
```

### Running pathway analysis
We can now use the [gage()](https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/gage) function to obtain the significantly perturbed pathways from our differential expression experiment. By default the [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) package performs this analysis while taking into account up and down regulation. Setting `same.dir=FALSE` will capture pathways perturbed without taking into account direction. This is generally not recommended for the GO groups as most genes within these gene sets are regulated in the same direction, however the same is not true for KEGG pathways and using this parameter may produce informative results in such cases.

```R
# Run enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(tumor_v_normal_DE.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(tumor_v_normal_DE.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(tumor_v_normal_DE.fc, gsets = go.cc.gs)

# covert the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# convert the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)
```

{% include question.html question="Which genes are in > 30% of significant pathways in the upregulated GO biological process results (q <= .05)" answer='Two genes are, RPS27A, UBA52. Here is an <a href="http://genviz.org/assets/pathway_analysis/exercise1/exercise1_pathway_analysis.R">Rscript</a> to get the correct answer.'%}

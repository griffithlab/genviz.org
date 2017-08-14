---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Pathway Analysis
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-01-01
---

In the previous section we examined differential expression of genes from the [E-GEOD-50760](https://www.ncbi.nlm.nih.gov/pubmed/25049118) data set. In this section we will use the [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) package to determine if there are any coordinated differential expression patterns in the data set we used for differential expression,  [E-GEOD-50760](https://www.ncbi.nlm.nih.gov/pubmed/25049118).

### What is gage?
generally applicable gene-set enrichment [(gage)](https://bioconductor.org/packages/release/bioc/html/gage.html) is a popular bioconductor package for performing gene-set and pathway analysis. The package works independent of sample sizes, experimental designs, assay platforms, and is applicable to both microarray and rnaseq data sets. In this section we will use [(gage)](https://bioconductor.org/packages/release/bioc/html/gage.html) and gene sets from the "Kyoto Encyclopedia of Genes and Genomes" [(KEGG)](http://www.kegg.jp/) and "Gene Ontology" [(GO)](http://www.geneontology.org/) databases to perform pathway analysis. First let's go ahead and install [(gage)](https://bioconductor.org/packages/release/bioc/html/gage.html) and load the differential expression results from the previous section.

```R
# install gage
source("https://bioconductor.org/biocLite.R")
biocLite("gage")

# load the differential expression results fro the previous section
load(url("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/deseq2Data_v1.RData"))
```

```R
################################################################################
############### annotate for GSEA ##############################################

library("AnnotationDbi")
library("org.Hs.eg.db")
require("gage")

# map the entrezid, gene name, and symbol to the deseq2 results
deseq2.res$symbol = mapIds(org.Hs.eg.db, keys=row.names(deseq2.res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
deseq2.res$entrez = mapIds(org.Hs.eg.db, keys=row.names(deseq2.res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
deseq2.res$name =   mapIds(org.Hs.eg.db, keys=row.names(deseq2.res), column="GENENAME", keytype="ENSEMBL", multiVals="first")

################################################################################
############## Compare to other papers #########################################

# chen et al.
chen <- read.delim("journal.chen.S5.txt")
deseq2.res$chen_et_al <- 0
deseq2.res[rownames(deseq2.res) %in% chen$Ensembl.Gene.ID,"chen_et_al"] <- 1

# output results
deseq2.res.df <- as.data.frame(deseq2.res)
deseq2.res.df <- deseq2.res.df[order(deseq2.res.df$padj),]
write.table(deseq2.res.df, file="DEseq2.res.tsv", sep="\t", row.names=T, quote=F)

################################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GSEA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
################################################################################

# grab the log fold changes for everything
deseq2.fc <- deseq2.res$log2FoldChange
names(deseq2.fc) <- deseq2.res$entrez
exp.fc <- deseq2.fc
out.suffix <- "deseq2"

# grab the log fold changes for just what is significant and unique in our data set compared to chen paper
deseq2.res.v2 <- deseq2.res[deseq2.res$pvalue < .05 & deseq2.res$chen_et_al == 0,]
deseq2.fc.v2 <- deseq2.res.v2$log2FoldChange
names(deseq2.fc.v2) <- deseq2.res.v2$entrez
exp.fc.v2 <- deseq2.fc.v2

################################################################################
################# Set up databases for enrichment analysis #####################
# set up kegg database
kg.hsa <- kegg.gsets("hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

################################################################################
################# Run the enrichment algorithms ################################

# Run enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(exp.fc, gsets = kegg.sigmet.gs, ref = NULL, samp = NULL)
fc.kegg.dise.p <- gage(exp.fc, gsets = kegg.dise.gs, ref = NULL, samp = NULL)
fc.go.bp.p <- gage(exp.fc, gsets = go.bp.gs, ref = NULL, samp = NULL)
fc.go.mf.p <- gage(exp.fc, gsets = go.mf.gs, ref = NULL, samp = NULL)
fc.go.cc.p <- gage(exp.fc, gsets = go.cc.gs, ref = NULL, samp = NULL)

# Run enrichment on genes unique to our dataset compared to chen paper
fc.kegg.sigmet.p.v2 <- gage(exp.fc.v2, gsets = kegg.sigmet.gs, ref = NULL, samp = NULL)
fc.kegg.dise.p.v2 <- gage(exp.fc.v2, gsets = kegg.dise.gs, ref = NULL, samp = NULL)
fc.go.bp.p.v2 <- gage(exp.fc.v2, gsets = go.bp.gs, ref = NULL, samp = NULL)
fc.go.mf.p.v2 <- gage(exp.fc.v2, gsets = go.mf.gs, ref = NULL, samp = NULL)
fc.go.cc.p.v2 <- gage(exp.fc.v2, gsets = go.cc.gs, ref = NULL, samp = NULL)

################################################################################
################ Output all results to tables ##################################

################ for all genes #################################################

# for kegg metabolism and signaling pathways
write.table(na.omit(fc.kegg.sigmet.p$greater), file="../gage/hcc_noncirrhotic.greater.fc.kegg.sigmet.p.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.kegg.sigmet.p$less), file="../gage/hcc_noncirrhotic.less.fc.kegg.sigmet.p.tsv", row.names=TRUE, sep="\t")

# for kegg disease pathways
write.table(na.omit(fc.kegg.dise.p$greater), file="../gage/hcc_noncirrhotic.greater.fc.kegg.dise.p.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.kegg.dise.p$less), file="../gage/hcc_noncirrhotic.less.fc.kegg.dise.p.tsv", row.names=TRUE, sep="\t")

# for gene ontologies: biological process
write.table(na.omit(fc.go.bp.p$greater), file="../gage/hcc_noncirrhotic.greater.fc.go.bp.p.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.go.bp.p$less), file="../gage/hcc_noncirrhotic.less.gc.go.bp.p.tsv", row.names=TRUE, sep="\t")

# for gene ontologies: molecular function
write.table(na.omit(fc.go.mf.p$greater), file="../gage/hcc_noncirrhotic.greater.fc.go.mf.p.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.go.mf.p$less), file="../gage/hcc_noncirrhotic.less.gc.go.mf.p.tsv", row.names=TRUE, sep="\t")

# for gene ontologies: cellular component
write.table(na.omit(fc.go.cc.p$greater), file="../gage/hcc_noncirrhotic.greater.fc.go.cc.p.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.go.cc.p$less), file="../gage/hcc_noncirrhotic.less.gc.go.cc.p.tsv", row.names=TRUE, sep="\t")

##################### Compared to chen paper ###################################

# for kegg metabolism and signaling pathways
write.table(na.omit(fc.kegg.sigmet.p.v2$greater), file="../gage/hcc_noncirrhotic.greater.fc.kegg.sigmet.p.v2.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.kegg.sigmet.p.v2$less), file="../gage/hcc_noncirrhotic.less.fc.kegg.sigmet.p.v2.tsv", row.names=TRUE, sep="\t")

# for kegg disease pathways
write.table(na.omit(fc.kegg.dise.p.v2$greater), file="../gage/hcc_noncirrhotic.greater.fc.kegg.dise.p.v2.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.kegg.dise.p.v2$less), file="../gage/hcc_noncirrhotic.less.fc.kegg.dise.p.v2.tsv", row.names=TRUE, sep="\t")

# for gene ontologies: biological process
write.table(na.omit(fc.go.bp.p.v2$greater), file="../gage/hcc_noncirrhotic.greater.fc.go.bp.p.v2.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.go.bp.p.v2$less), file="../gage/hcc_noncirrhotic.less.gc.go.bp.p.v2.tsv", row.names=TRUE, sep="\t")

# for gene ontologies: molecular function
write.table(na.omit(fc.go.mf.p.v2$greater), file="../gage/hcc_noncirrhotic.greater.fc.go.mf.p.v2.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.go.mf.p.v2$less), file="../gage/hcc_noncirrhotic.less.gc.go.mf.p.v2.tsv", row.names=TRUE, sep="\t")

# for gene ontologies: cellular component
write.table(na.omit(fc.go.cc.p.v2$greater), file="../gage/hcc_noncirrhotic.greater.fc.go.cc.p.v2.tsv", row.names=TRUE, sep="\t")
write.table(na.omit(fc.go.cc.p.v2$less), file="../gage/hcc_noncirrhotic.less.gc.go.cc.p.v2.tsv", row.names=TRUE, sep="\t")

################################################################################
########### Determine what is driving significant pathways #####################

# grab data and make sure pathways significant
keggGreater <- as.data.frame(fc.kegg.sigmet.p$greater)
keggLess <- as.data.frame(fc.kegg.sigmet.p$less)
goGreater <- as.data.frame(fc.go.bp.p$greater)
goLess <- as.data.frame(fc.go.bp.p$less)

keggGreater <- na.omit(keggGreater[keggGreater <= .05,])
keggLess <- na.omit(keggLess[keggLess <= .05,])
goGreater <- na.omit(goGreater[goGreater <= .05,])
goLess <- na.omit(goLess[goLess <= .05,])

# get genes associated with pathways
a <- function(x, y){
    xID <- rownames(x)
    xEntrez <- kegg.sigmet.gs[which(names(y) %in% xID)]
    xEntrez <- unlist(xEntrez)
    xEntrez <- plyr::count(xEntrez)
    xEntrez$symbol <- mapIds(org.Hs.eg.db, keys=as.character(xEntrez$x), column="SYMBOL", keytype="ENTREZID", multiVals="first")
    return(xEntrez)
}

keggGreaterCount <- a(keggGreater, kegg.sigmet.gs)
keggLessCount <- a(keggLess, kegg.sigmet.gs)
goGreaterCount <- a(goGreater, go.bp.gs)
goLessCount <- a(goLess, go.bp.gs)
```


### KEGG
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### GO
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

### GAGE
Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

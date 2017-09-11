---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Variant annotation and interpretation
categories:
    - Module 5
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-01-01
---

When variants are identified in the genome (or transcriptome) some kind of annotation and need for interpretation invariably follows. There are many, many tools for annotation and interpretation in different contexts and for different purposes. In this section we explore just a few of these many options. First we will learn to use Ensembl's Variant Effect Predictor (VEP), a popular and widely used variant transcript annotator. VEP has many functions, but it is first used to annotate variants in the context of set of known transcripts. The other resources we will use, ClinVar and CIViC attempt to summarize evidence for the clinical relevance of variants in inherited human diseases and cancer respectively.

Some here are some examples of variant annotation and interpretation contexts:
* Population frequency/recurrence (is the variant common, rare, rare?)
* Transcripts (does the variant occur within a transcribed region of a gene? Does is affect the predicted translation of that transcript?)
* Function (is the variant likely to disrupt the normal function of a gene?). There are many, many approaches to this.
  * Conservation of the affected region
  * Predicted biochemical significance of amino acid alterations
  * Occurence in know functional domains (e.g. the binding pocket of a kinase)
  * Hot spots of variantion (some patterns can suggest gain-of-function)
    * 2D hotspots
    * 3D hotspots
  * Patterns that suggest loss of function
  * Actual experimental evidence for the specific variant or one very similar

* What other approaches can you think of?

[Module 5 Lecture](https://github.com/griffithlab/gen-viz-lectures/raw/master/GenViz_Module5_Lecture.pdf)

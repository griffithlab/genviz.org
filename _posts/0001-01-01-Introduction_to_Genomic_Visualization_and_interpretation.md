---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to Genomic Data Visualization and Interpretation
categories:
    - Module 1
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-01-01
---

Welcome to the Genomic Visualization and Interpretations tutorial series. The goal of this tutorial is to provide students with the skills necessary to visualize and interpret 'omic' data in various forms. To accomplish this task we will be exploring a number of tools and resources available through the R programming language and via the web providing background where appropriate. We recommend using the navigation links found at the bottom of all course topic pages however students can search for individual topics via the search tab. Exercises to reinforce concepts are interspersed throughout each course topic.

***
## Central Dogma

{% include figure.html image="/assets/central_dogma.png" position="right" %}

> *"The Central Dogma. This states that once 'information' has passed into protein it cannot get out again. In more detail, the transfer of information from nucleic acid to nucleic acid, or from nucleic acid to protein may be possible, but transfer from protein to protein, or from protein to nucleic acid is impossible. Information means here the precise determination of sequence, either of bases in the nucleic acid or of amino acid residues in the protein."*
> --Francis Crick 1956

As quoted above the central dogma describes the flow of biological sequence information encoded in deoxyribonucleic acid (DNA), ribonucleic acid (RNA) and protein. In essence it states that information cannot flow backward from a protein, that is to say transfer of information from a protein to RNA, DNA, or replicated by another protein is impossible. The general flow of information is through the copying of DNA through DNA replication, the creation of RNA from transcription of DNA,  and the creation of proteins via translation of RNA. This process occurs in most cells and are considered mechanisms of general transfer. Special transfers do exist in the form of RNA replication and reverse transcription however these processes are mostly limited to viruses.

Delving in a bit deeper, DNA takes the form of a double stranded helix comprised of pairs of complementary bases. Adenine (A) and Guanine (G) are classified as purines and complement the pyrimidines Thymine (T) and Cytosine (C) respectively. During transcription DNA is transcribed into a single stranded mRNA molecule in which T bases are converted to uracil (U) and RNA editing such as the removal of introns is performed. This results in a mature mRNA structure composed of a 5' UTR, Coding sequence (CDS), 3' UTR, and a polyadenylated tail. The mature mRNA is then translated into a peptide sequence beginning at the AUG start codon within the CDS and folded into a protein.

{% include question.html question="What is the mRNA sequence of: ATGTTTACTGCTGATGGCCGCTGA?" answer="AUGUUUACUGCUUGAUGGCCGCUGA"%}
{% include question.html question="What is the translated peptide sequence in the previous question?" answer="Met-Phe-Thr-Ala-Asp-Gly-Arg-Opal"%}

***

## Omic technologies and Data

The advent of rapid and cheap massively parallel sequencing has dramatically increased the availability of genome, transcriptome, and epigenome data. Further this has given rise to standardized workflows and file formats. Extracting biologically meaningful data and interpreting the results remains challenging.

### genomics



### transcriptomics

### epigenomics

 ***

## Reference Files

Over this course we will be working with a number of refernce files to aid in analysis and interpretation of our data. All of these files are in standardized formats and are freely available. A brief description of each reference file is given below.

- [FASTA](http://genetics.bwh.harvard.edu/pph/FASTA.html){:target="_blank"}: Text files which contain nucleotide or peptide sequences.
- [GTF](http://www.ensembl.org/info/website/upload/gff.html){:target="_blank"}: (Gene Transfer Format) Tab delimited files which hold information regarding gene structure.
- [BED](http://www.ensembl.org/info/website/upload/bed.html){:target="_blank"}: (Browser Extensible Data) Tab delimited files hold feature information commonly used to diplay data on an annotation track. Cordinates within these files are 0-based.
- [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf){:target="_blank"} (Variant Call Format) Text file used to store observed sequence variations.
- [VEP](http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#output){:target="_blank"} (Variant Effect Predictor) Annotation file used to provide additional information for sequence variations. Originates from the ensembl VEP algortihm.
- [MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification){:target="_blank"} (Mutation Annotation Format) Annotation file used to provide additional information for sequence variations. Widely used in the cancer genome atlas.
- [SAM/BAM/CRAM](https://samtools.github.io/hts-specs/SAMv1.pdf) Sequence alignment map format and it's compressed equivalents binary and compressed alignment map are files for storing aligned sequencing data. BAM files are just binary files of the SAM file. CRAM uses the reference sequence to more effeciently compress the information in a SAM file. These file types are commonly viewed and manipulated with [samtools](https://github.com/samtools/samtools).

{% include question.html question="What character designates a header in a fasta file?" answer="\">\""%}
{% include question.html question="What is the difference between 0-based and 1-based coordinates?" answer="0-based coordiantes number between nucleotides, 1-based coordinates number nucleotides directly."%}

***

## Genomic annotation resources, browsers, etc.

Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

***

## Introduction to demonstration data settings

Lorem ipsum dolor sit amet, munere intellegat cu mel. Ea sint summo exerci mei. Autem tritani scaevola mei ea, sonet oporteat vel cu. Duo cu erat libris vulputate. Cum possim copiosae facilisi ea, partiendo tincidunt voluptatibus ne est, vix ea justo animal.

Cum quem justo urbanitas no, mei inermis alienum indoctum ei. Cu assum ludus soluta per. Sea at idque perpetua, ex fabulas hendrerit adversarium per, sit impedit recteque necessitatibus an. Quo fabulas feugait scriptorem et.

***

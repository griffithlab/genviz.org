---
feature_text: |
  ## Genomic Visualization and Interpretations
title: Introduction to BioMart
categories:
    - Module 1
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-06-01
---

[Ensembl BioMart](http://www.ensembl.org/biomart/martview) is a powerful web tool (with API) for performing complex querying and filtering of the various Ensembl databases (Ensembl Genes, Mouse Strains, Ensembl Variation, and Ensembl Regulation). It is often used for ID mapping and feature extraction. Almost any data that is viewable in the Ensembl genome browser can be accessed systematically from BioMart. 

### Browsing Ensembl BioMart
To browse [Ensembl BioMart](http://www.ensembl.org/biomart/martview) you can simply navigate to the main page directly, by google or from any Ensembl page. 

{% include figure.html image="/assets/BioMart/BioMart_start.png" width="800" %}

You will first need to CHOOSE DATABASE which can be any of Ensembl Genes, Mouse Strains, Ensembl Variation, or Ensembl Regulation. Note, the database version is included in the name (e.g., Ensembl Genes 90). The Ensembl Genes dataset is probably the most commonly desired choice, depending on your purpose. Next, you will CHOOSE DATASET, which for Genes means selecting a species (e.g., Human genes (GRCh38.p10)). For illustration, we will walk through some examples using the **Ensembl Genes 90** database and **Human genes** dataset. At this point you have the option to apply filters and desired attributes from left panel. Note that two attributes (Gene stable ID, and Transcript stable ID) have been pre-selected for us.

{% include figure.html image="/assets/BioMart/BioMart_prefilter.png" width="800" %}




### Using older versions of Ensembl BioMart
By default, [Ensembl BioMart](http://www.ensembl.org/biomart/martview) only presents data for the latest, most current version of Ensembl. Older versions can be accessed by navigating to an [Ensembl Archive Site](http://www.ensembl.org/Help/ArchiveList) (linked from the bottom right of every Ensembl page) and then following the BioMart link (top left of every Ensembl Archive page). For example, the last version of Ensembl for the human GRCh37 (hg19) build was v75 (February 2014). The Archive EnsEMBL release 75 was available at [http://feb2014.archive.ensembl.org](http://feb2014.archive.ensembl.org) and the corresponding BioMart at [http://feb2014.archive.ensembl.org/biomart/martview](http://feb2014.archive.ensembl.org/biomart/martview).


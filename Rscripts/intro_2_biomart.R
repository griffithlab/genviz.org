# Load the BioMart library
library("biomaRt")

# Create a mart object from ensembl Biomart
mart <- useMart("ensembl", "hsapiens_gene_ensembl")

# Output all attributes to see which are available and how they are named
attributes <- listAttributes(mart)
attributes

# Output all filters to see which are available and how they are named
filters <- listFilters(mart)
filters

# Get Ensembl gene records and relevant attributes for a list of HGNC symbols
hgnc_symbols <- c("NEUROD2", "PPP1R1B", "STARD3", "TCAP", "PNMT", "PGAP3", "ERBB2", "MIR4728", "MIEN1", "GRB7", "IKZF3")
annotations_ENSG=getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","external_gene_name","hgnc_symbol"), filter="hgnc_symbol", values=hgnc_symbols, mart=mart)

################################################################################
################ Exercise 1 ####################################################

#!!!!!!!!!! "fc.go.bp.p.up" and "go.bp.gs" should already be in your R environment, see tutorial above !!!!!!!!!#

# get genes associated with pathways
a <- function(sigPathways, geneSetDB){
    
    # grab the names of all significant pathways
    sigPathwaysID <- rownames(sigPathways)
    
    # subset the geneSet list to only those in the significant pathways
    sigPathwaysGenes <- geneSetDB[which(names(geneSetDB) %in% sigPathwaysID)]
    numberSigPathways <- length(sigPathwaysGenes)
    sigPathwaysGenes <- unlist(sigPathwaysGenes)
    
    # count the number of times a gene occurs in these significant pathways
    sigPathwaysTable <- plyr::count(sigPathwaysGenes)
    
    # annotate these final results with the gene symbol and some extra information
    sigPathwaysTable$symbol <- mapIds(org.Hs.eg.db, keys=as.character(sigPathwaysTable$x), column="SYMBOL", keytype="ENTREZID", multiVals="first")
    sigPathwaysTable$proportion <- sigPathwaysTable$freq/numberSigPathways
    
    # sort and return the results
    sigPathwaysTable <- sigPathwaysTable[order(-sigPathwaysTable$proportion),]
    
    return(sigPathwaysTable)
}

# run the function defined above
geneTable <- a(na.omit(fc.go.bp.p.up[fc.go.bp.p.up$q.val <= .05,]), go.bp.gs)
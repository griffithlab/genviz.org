# load the ouput of samtools depth into R
seqData <- read.delim("~/Desktop/ALL1_CaptureDepth.tsv", header=F)
colnames(seqData) <- c("chr", "position", "samp1", "samp2", "samp3", "samp4")
seqData <- seqData[,c(3:6)]

# Count the occurrences of each coverage value
seqCovList <- apply(seqData, 2, count)

# rename the columns for each dataframe
renameCol <- function(x){
    colnames(x) <- c("coverage", "freq")
    return(x)
}
seqCovList <- lapply(seqCovList, renameCol)

# create framework data frame with entries for the min to the max coverage
maximum <- max(unlist(lapply(seqCovList, function(x) max(x$coverage))))
minimum <- min(unlist(lapply(seqCovList, function(x) min(x$coverage))))
covFramework <- data.frame("coverage"=minimum:maximum)

# Merge the framework data frame with the coverage
# seqCovList <- lapply(seqCovList, function(x, y) merge(x, y, by="coverage", all=TRUE), covFramework)
seqCovList <- lapply(seqCovList, merge, covFramework, by="coverage", all=TRUE)

# merge all data frames together
seqCovDataframe <- Reduce(function(...) merge(..., by="coverage", all=T), seqCovList)

# set all NA values to 0
seqCovDataframe[is.na(seqCovDataframe)] <- 0

# set the rownames, remove the extra column, and convert to a matrix
rownames(seqCovDataframe) <- seqCovDataframe$coverage
seqCovDataframe$coverage <- NULL
seqCovMatrix <- as.matrix(seqCovDataframe)

# rename columns
colnames(seqCovMatrix) <- c("sample1", "sample2", "sample3", "sample4")

# run covbars
covBars(seqCovMatrix)
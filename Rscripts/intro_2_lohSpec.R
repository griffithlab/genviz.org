################################################################################
############################ Tutorial ##########################################
################################################################################
# read in the varscan data 
lohData <- read.delim("~/Desktop/HCC1395.varscan.tsv", header=FALSE)

# grab only those columns which are required and name them
lohData <- lohData[,c("V1", "V2", "V7", "V11")]
colnames(lohData) <- c("chromosome", "position", "n_vaf", "t_vaf")

# convert the normal and tumor vaf columns to fractions
lohData$n_vaf <- as.numeric(gsub("%", "", lohData$n_vaf))/100
lohData$t_vaf <- as.numeric(gsub("%", "", lohData$t_vaf))/100

# add a sample column
lohData$sample <- "HCC1395"

# limit to just the genome of interest
lohData <- lohData[grepl("^\\d|X|Y", lohData$chromosome),]

# create an inital plot
lohSpec(lohData)

# Obtain variants with a VAF greater than 0.4 and below 0.6.
lohData <- lohData[lohData$n_vaf > 0.4 & lohData$n_vaf < 0.6,]

# run lohSpec
lohSpec(x=lohData)

################################################################################
############################ Exercises #########################################
################################################################################
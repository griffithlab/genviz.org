
# Read in data
data <- read.table(file="~/Downloads/testset_gcrma.txt", header=TRUE, sep="\t", as.is=c(1:3))

# Extract just the expression data
expression_matrix <- data[,4:length(colnames(data))]

# Convert expression data from log2 to raw intensity
raw_expression_matrix <- 2^expression_matrix

# Create a custom function to return % values >= cutoff
my_pe_fun <- function(x, cutoff){
  perc_expressed <- (length(which(x >= cutoff))/length(x))*100
  return(perc_expressed)
}

# Create a custom function to return coefficient of variation value
my_cov_fun <- function(x){
  coeff_variation <- sd(x)/mean(x)
  return(coeff_variation)
}

# Apply the above functions to all probes in the raw intensity matrix
percent_expressed <- apply(raw_expression_matrix, 1, my_pe_fun, 100)
cov <- apply(raw_expression_matrix, 1, my_cov_fun)

#Create a new dataframe that includes the percent expressed and cov values that we will filter on
data_plus=cbind(data,percent_expressed,cov)

# Extract subset of data matrix that pass both conditions
filtdata = data_plus[data_plus$percent_expressed>=20 & data_plus$cov>=0.7,]

#Calculate the number of probesets meeting the filter criteria
length(filtdata[,"probes"])
dim(filtdata)

# exercise 1
# create a custom genome
chr10 <- data.frame("chromosome"="chr10", "start"=0, "end"=135534747)

# create custom layer to highlight q23.1
library(ggplot2)
layer1 <- geom_vline(xintercept=c(89500000, 92900000), colour="chartreuse", linetype=2, size=1)

# create the plot
lohSpec(lohData[lohData$chromosome==10,], y=chr10, plotLayer = layer1)
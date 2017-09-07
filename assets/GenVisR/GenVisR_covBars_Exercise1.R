################################################################################
################# Exercise 1 ###################################################
################################################################################

# load ggplot2
library(ggplot2)

############### Make additional layers to pass to covBars ######################

# alter the y-axis facet labels
colnames(seqCovMatrix2) <- c("Skin_d42_I\n(skin normal)", "SL_d3072_I\n(sorted lymphs)", "SB_d3072_A\n(relapse 2 sorted blasts)", "M_d3068_A\n(relapse 2 bluk marrow)")

# remove the grey facet boxes
layer1 <- theme(strip.background=element_rect(fill="white"))

# position the facet text on the right side
layer2 <- facet_grid(sample ~ ., switch="y")
layer3 <- theme(strip.text.y=element_text(angle=180))

# add a plot title
layer4 <- ggtitle("Custom Capture Data")

# change the x-axis title
layer5 <- xlab("Cumulative Coverageof Region\nTargeted by and Successfully\nTiled by Nimblegen")

# change the legend
rainbowColorRamp <- rainbow(1200)[1:1050]
layer6 <- scale_fill_gradientn(name="Seq.\nDepth", colours=rainbowColorRamp, labels=c("0", "300", "600", "900", "1200+"))

# run covBars with all of these layers
covBars(seqCovMatrix2, plotLayer = list(layer1, layer2, layer3, layer4, layer5, layer6))
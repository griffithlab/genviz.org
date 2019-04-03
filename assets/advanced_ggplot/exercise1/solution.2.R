# save finalGrob as exercise 1
exercise1 <- finalGrob

# change the legend border to purple for all 3 legends
finalGrob$grobs[[1]]$grobs[[15]]$grobs[[1]]$grobs[[1]]$gp$col <- "slateblue3"
finalGrob$grobs[[2]]$grobs[[15]]$grobs[[1]]$grobs[[1]]$gp$col <- "slateblue3"
finalGrob$grobs[[3]]$grobs[[15]]$grobs[[1]]$grobs[[1]]$gp$col <- "slateblue3"

# plot the result
grid.draw(finalGrob)
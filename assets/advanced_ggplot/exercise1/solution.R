# define layout
layout <- rbind(c(1, 1, 1, 1, 1),
                c(1, 1, 1, 2, 2),
                c(NA, NA, NA, NA, NA),
                c(NA, 4, 4, 3, NA))

# make plot
grid.arrange(grob1, grob2, grob3, grob4, layout_matrix=layout)

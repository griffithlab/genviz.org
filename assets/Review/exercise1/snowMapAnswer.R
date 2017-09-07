# Install and load packages
#install.packages("cholera")
library(cholera)
library(ggplot2)

# load data frames from the cholera package
data(roads)
data(pumps)
data(fatalities.address)
data(pump.case)

# transform pump.case to a df
pump.case.df <- stack(pump.case)
colnames(pump.case.df) <- c("anchor.case", "pump")

# remove the p designation
pump.case.df$pump <- gsub("p", "", pump.case.df$pump)
pump.case.df <- merge(fatalities.address, pump.case.df, by="anchor.case")
pump.case.df <- merge(pump.case.df, pumps[,c("id", "street")], by.x="pump", by.y="id")

# ggplot2 call
ggplot() + geom_line(data=roads, aes(x = x, y = y, group=street), colour="white") +
    geom_point(data=pumps, aes(x=x, y=y, shape="pump"), colour="red", size=3) +
    scale_shape_manual("Pump Location", values=17) +
    geom_point(data=pump.case.df, aes(x=x, y=y, colour=street)) +
    scale_colour_manual(values=c("#a2ffe2", "#dda3f1", "#85bd70", "#f3daff",
                                 "#d5ed93", "#b6a5c3", "#26eeeb", "#f49f7d",
                                 "#b6bdae", "#acb069")) +
    geom_density2d(data=pump.case.df, aes(x=x, y=y), colour="grey70") +
    guides(colour=guide_legend(title="Deaths By Nearest Pump")) +
    coord_fixed() + theme_void() +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.background=element_rect(fill="grey40")) + ggtitle("1854 Cholera Epedemic")

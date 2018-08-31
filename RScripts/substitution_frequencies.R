# Script called from python lib called VarPlotLib
# Used to plot base substitution frequencies by groups
# Author: Yrj?? Koski, 2018

library(ggplot2)
library(grid)
library(RColorBrewer)

inputArgs <- commandArgs(trailingOnly = TRUE)
input <- inputArgs[1]

mutations <- read.csv(input, header=TRUE, sep=",", check.names = FALSE)
n <- length(levels(mutations$Mutation))

getPalette = colorRampPalette(brewer.pal(11,"Spectral"))

g <- ggplot(mutations) + aes(x = Group, fill = factor(Mutation), weight = Amount) + geom_bar(position = "fill") +
  scale_fill_manual(values = getPalette(n)) + theme(legend.title=element_blank()) + 
  scale_y_continuous(expand = c(0,0))

png(file="Mutation_freqs.png",width = 200, height = 600)
  grid.draw(g)
dev.off()


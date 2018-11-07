# Script called from python lib called VarPlotLib
# Used to plot variants by genes in multiple samples
# Author: Yrjo Koski, 2018

library(ggplot2)
library(reshape2)
library(grid)
library(RColorBrewer)
library(dplyr)
library(gtable)
library(gridExtra)

'%ni%' <- Negate('%in%')

rename_samples <- function(group.name, matrix) {
  c <- 1
  for (i in 1:nrow(matrix)){
    if (as.character(matrix[i,"Group"]) == group.name) {
      matrix[i,"Sample"] <- paste(group.name,'_',as.character(c),sep = '')
      c <- c + 1
    }
  }
  return(matrix)
}

inputArgs <- commandArgs(trailingOnly = TRUE)
gene.matrix <- inputArgs[1]
group.table <- inputArgs[2]
var.amount <- inputArgs[3]
custom.height <- inputArgs[4]
custom.width <- inputArgs[5]
rename.samples <- inputArgs[6]
removeNA <- inputArgs[7]


var.mat <- read.csv(gene.matrix, sep=',', header = TRUE)
# Read amount of variants per sample from var.amount
n.vars <- read.csv(var.amount, sep=',', header = TRUE)
gene.groups <- read.table(group.table, sep='\t', header = FALSE)
if (length(gene.groups[1,]) == 1) {
  gg = FALSE
} else {
  gg = TRUE
}

var.mat[var.mat=="-"] <- NA
var.mat$Sample <- as.character(var.mat$Sample)
n.vars <- cbind(n.vars, Group=var.mat$Group)

# Remove rows with only NA-values
# Do the same for variant amounts
if (removeNA == TRUE) { 
  var.mat <- var.mat[rowSums(is.na(var.mat)) != ncol(var.mat)-2, ]
  n.vars <- n.vars[n.vars$Sample %in% var.mat$Sample, ]
}

if (rename.samples == TRUE) {
  for (i in levels(var.mat$Group)) {var.mat <- rename_samples(i,var.mat)}
}

# Sort df by amount of variant in each row
var.mat.sorted <- var.mat[order(rowSums(is.na(var.mat))), ]
var.mat.sorted <- var.mat.sorted[,order(colSums(is.na(var.mat.sorted)), decreasing = TRUE)]

# Set height and width of the plot. If custom values are given use them, else approximate them
if (custom.height != FALSE) {h <- as.integer(custom.height)} else {h <- 120 + (ncol(var.mat.sorted) - 2) * 30}
if (custom.width != FALSE) {w <- as.integer(custom.width)} else {w <- 80 + (nrow(var.mat.sorted)) * 30}
#str(var.mat.sorted)

var.mat.melt <- melt(var.mat.sorted, id.vars = c('Sample','Group'))
var.mat.melt$value <- as.character(var.mat.melt$value)

var.types <- c("m","f","n","s","d",NA)
var.mat.melt$value[var.mat.melt$value %ni% var.types] <- "mul"
var.mat.melt$value <- factor(var.mat.melt$value, levels = c("m","f","n","s","d","mul"))

# Plot variant matrix
xlims <- c(0,length(var.mat$Sample))
g <- ggplot(var.mat.melt, aes(Sample, variable)) + geom_tile(aes(fill = value), size = 0.2, color = "white")
      if (gg == TRUE) {g <- g + facet_grid(gene.groups$V2[variable] ~ ., scales = "free", space = "free", switch = "y")
      } else {
          g <- g + facet_grid( . ~ Group, scales = "free", space = "free", switch = "y")
      }

g <- g + scale_fill_manual(name="Variant",
                           breaks=c("m","f","n","s","d","mul"),
                           labels=c("Missense","Frameshift","Nonsense","Splicing_site","Nonframeshift_deletion","Multiple"),
                           values = c("#ABDDA4", "#D53E4F", "#FEE08B", "#3288BD", "#F46D43", "#5E4FA2"),#, "#8C510A"), 
                           drop = FALSE) +
  scale_x_discrete(position = "top")

if (length(is.na(var.mat.melt$Group)) == length(var.mat.melt$Group)) 
  {g <- g + theme(strip.background = element_blank(),
                  strip.text.x = element_blank(),
                  axis.text.x = element_text(angle = -90, hjust = 1),
                  plot.title = element_text(hjust = 0.5),
                  panel.grid.major = element_blank()) + labs(y='')
  } else {
    g <- g + theme(axis.text.x = element_text(angle = -90, hjust = 1),
              plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank()) +labs(y='')
    }

  #coord_equal()

#n.vars$Amount_of_variants <- factor(n.vars$Amount_of_variants)

# Plot bar_graph
g2 <- ggplot(n.vars, aes(x = Sample)) + geom_bar(aes(y = Amount_of_variants), stat="identity", size = 0.2, color ="white") +
  facet_grid( ~ Group, scales = "free", space = "free_x") +
  theme(panel.grid.major = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        aspect.ratio = 6) +
  ylab("Amount of variants")
  guides(fill=FALSE)
 
p1 <- ggplotGrob(g)
p2 <- ggplotGrob(g2)
colnames(p1) <- paste0(seq_len(ncol(p1)))
colnames(p2) <- paste0(seq_len(ncol(p2)))

# grid.draw(gtable_combine(p1, p2, along=2))

png(file="Variant_matrix.png",width = w, height = h)
grid.draw(gtable_combine(p1, p2, along=2))
dev.off()

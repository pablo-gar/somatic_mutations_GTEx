# Creates a bar plot of the number of tissues to analyzed and to be analyzed
# Rscript tissues.R ~/results/FraserLab/somaticMutationsProject/20180220/tissues.pdf

library("rjson")
library("ggplot2")
library("reshape")
source("../R/ggthemes.R")

outFile <- commandArgs(T)[1]

dir.create(dirname(outFile))
if(!dir.exists(dirname(outFile)))
    stop("Could not create out folder")



colors <- c("#D9AC91", "#91BED9")

config <- fromJSON(file = "../config.json")

tissueFile <- config$auxiliaryFiles$tissues

tissues <- read.table(tissueFile, sep = "\t", header = T, stringsAsFactors = F)
tissues$NotProcessed <- tissues$Total.Samples - tissues$Processed
tissues <- tissues[,-2]
tissues <- melt.data.frame(tissues, id.vars = "Tissue")
tissues$variable <- factor(tissues$variable, levels = sort(unique(as.character(tissues$variable))), ordered = T)


#plot
p <- ggplot(tissues, aes(x = Tissue, y = value, fill = variable)) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    ylab("Samples") + xlab("Tissue") +
    theme_grid_y() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top", plot.margin= unit(c(0.5, 0.5, 0.5, 2), units = "cm") )
    
ggsave (outFile, p, width = 12, height = 4)
    

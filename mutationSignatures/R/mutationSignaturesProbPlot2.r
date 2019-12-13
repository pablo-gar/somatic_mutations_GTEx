## Given a mutation signature table (table W from NMF) 
## Plots the mutation probability distribution for each mutation signature
## present in the table

## Author: Pablo E. Garcia-Nieto
## Date: 12/4/1991

library("ggplot2")
library("reshape")
source("../../R/ggthemes.R")

##----------------------------------
## ARGS and GLOBAL
args <- commandArgs(T)

signatureFile <- args[1]
outFilePlot <- args[2]

signatureFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cancerSignaturesStandard.txt"
outFilePlot <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/consensusMutationSignatures/cancerSignatures2.pdf"

MUTATION_KEY <- c(`00` = "T>A", `01` = "T>T", `02` = "T>G", `03` = "T>C",
                  `10` = "C>A", `11` = "C>T", `12` = "C>G", `13` = "C>C")

barColors <- c("#3BC9F3", "#2B2E34", "#FC3218", "#CAD0CE", "#9CD169", "#F1CAC9")

HEIGHT = 2.2
WIDTH = 4
NCOL_PLOT = 6


##----------------------------------
## MAIN

signatureTable <- read.table(signatureFile, sep = "\t", stringsAsFactors = F, header = F, row.names = 1)
colnames(signatureTable) <-paste("Signature", 1:ncol(signatureTable))
signatureTable <- signatureTable[!grepl("01", rownames(signatureTable)) & !grepl("13", rownames(signatureTable)) , ]

# Calculate frequencies
signatureTable <- as.data.frame(t(t(signatureTable) / colSums(signatureTable)))
#signatureTable <- as.data.frame(consensusSignatures)

# Append context and mutation type
mutationCode <- do.call(rbind, strsplit(rownames(signatureTable), "\\."))
signatureTable$mut <- MUTATION_KEY[mutationCode[,1]]
signatureTable$context <- mutationCode[,2]
signatureTable <- signatureTable[signatureTable$context != "C>C" & signatureTable$context != "T>T", ]

# Re arrange for plotting
toPlot <- melt(signatureTable, id.vars = c("mut", "context"))
toPlot$x <- factor(paste0(toPlot$mut, ".", toPlot$context), 
                   levels = sort(unique(paste0(toPlot$mut, ".", toPlot$context))), 
                   ordered = T)

#toPlot$context <- factor(as.character(toPlot$context), levels = unique(as.character(toPlot$context)), ordered = T)

p <- ggplot(toPlot, aes(x = context, fill = mut, y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(variable~mut, ncol = NCOL_PLOT, scales = "free") + 
    theme_grid_y() +
    scale_fill_manual(values = barColors) + 
    ylab("Mutation probability") + 
    #coord_cartesian(ylim = c(0,0.2)) + 
    theme(axis.text.x = element_text(angle=45, hjust=1))

nSignatures <- length(unique(toPlot$variable))

nsign <- length(unique(toPlot$variable))
NROW_PLOT <- nsign
ggsave(outFilePlot, p, height = HEIGHT * NROW_PLOT, width = WIDTH * NCOL_PLOT, limitsize=F)
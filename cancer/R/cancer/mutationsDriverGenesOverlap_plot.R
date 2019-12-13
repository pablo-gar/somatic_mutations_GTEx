# Calculates the signifcance of the number of mutations in driver genes based
# on a comparison to a set of genes with similar expression patterns


# Usage
# Rscript mutationsDriverGenesOverlap_plot.R outPlot.pdf permutationTable1 [permutationTable2] [...]

library("tidyr")
library("dplyr")
library("ggplot2")
source("../../../R/misc.R")
source("../../../R/plots.R")
source("../../../R/ggthemes.R")

cmdArgs <- commandArgs(T)
outPlot <- cmdArgs[1]
permutationTables <- normalizePath(cmdArgs[2:length(cmdArgs)])

tissues <- basename(dirname(dirname(permutationTables)))

HEIGHT_POINT = 0.4
WIDTH_PLOT = 7


# To test
#tissues <- list.dirs("/scratch/users/paedugar/somaticMutationsProject/cancer/driverEnrichment", recursive = F)
#permutationTables <- file.path(tissues, "n6_0.0_0.5", "permutationsWilcoxResults.txt")
#outPlot <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverEnrichmentPlots/allDrivers/meanFC_overPermutation.pdf"
#tissues <- basename(tissues)
#tissues <- tissues[file.exists(permutationTables)]
#permutationTables <- permutationTables[file.exists(permutationTables)]

#----------------------------------#
# METHODS
#----------------------------------#

readPermutationTables <- function(tissues, tables) {
    
    results <- list()
    for (i in 1:length(tissues)) {
        x <- read.table(tables[i], sep = "\t", stringsAsFactors = F, header = T)
        x$tissue <- tissues[i]
        results[[i]] <- x
    }
    
    results <- do.call(rbind, results) 
    return(results)
}

readTestResults <- function(tissues, tables) {
    
    results <- list()
    for (i in 1:length(tissues)) {
        x <- read.table(tables[i], sep = "\t", stringsAsFactors = F, header = T)
        x$tissue <- tissues[i]
        results[[i]] <- x
    }
    
    results <- do.call(rbind, results) 
    return(results)
}
    
    

#----------------------------------#
# MAIN
#----------------------------------#

#----------#
# Doing plots for mean mutation fold change values of driver genes and permuted genes (same expression)

permResults <- readPermutationTables(tissues, permutationTables)
permResultsGathered <- gather(permResults, "type", "log2Fold_overPermutation", realMeanDiff, randomMeanDiff)

# Calculating mean differences
toPlot <-
    permResultsGathered %>%
    group_by(tissue, type) %>%
    summarise(mean = mean(log2Fold_overPermutation), sd = sd(log2Fold_overPermutation)) %>%
    mutate(maxError = mean + sd,
           minError = mean - sd
           ) %>%
    ungroup()

# Getting wilcox pvalues
toPlotSignifLabels <-
    permResults %>%
    group_by(tissue) %>%
    summarise(wilcox.pvalue = wilcox.test(realMeanDiff, randomMeanDiff)[["p.value"]]) %>%
    mutate(wilcox.pvalueBon = p.adjust(wilcox.pvalue, method = "bonf")) %>%
    mutate(labelSignif = labelPvalues(wilcox.pvalueBon)) %>%
    ungroup()

# ordering by diference of means random vs real
meanDiff <- tapply(toPlot$mean, toPlot$tissue, function(x) x[1] - x[2])
toPlot$tissue <- factor(toPlot$tissue, levels = names(sort(meanDiff)), ordered = T)

p <- pointRange(as.data.frame(toPlot), x = "tissue", y = "mean", colour = "type", errorBarYmax = "maxError", errorBarYmin = "minError", pSize = 0.5) + 
         geom_text(aes(label = labelSignif), data = toPlotSignifLabels, y = max(toPlot$maxError)) +
         ylab("Mean FC (drivers/random)") +
         coord_flip() +
         scale_colour_manual(values = c("grey70", "grey20")) +
         theme_grid_x()
     
     
ggsave(outPlot, p, width = WIDTH_PLOT, height = length(unique(toPlot$tissue)) * HEIGHT_POINT)

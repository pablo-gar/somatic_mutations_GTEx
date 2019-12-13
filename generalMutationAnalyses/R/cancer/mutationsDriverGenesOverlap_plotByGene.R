# Calculates the signifcance of the number of mutations in driver genes based
# on a comparison to a set of genes with similar expression patterns

# Usage
# Rscript mutationsDriverGenesOverlap_plotByGene.R outPlot.pdf outPlot2.pdf permutationTable1 [permutationTable2] [...]

library("tidyr")
library("dplyr")
library("ggplot2")
library("gplots")
source("../../../R/misc.R")
source("../../../R/plots.R")
source("../../../R/ggthemes.R")
source("../../../R/geneTools.R")

cmdArgs <- commandArgs(T)
heatmapOut <- parseArg(cmdArgs[1], ",")
histOut <- parseArg(cmdArgs[2], ",")
permutationTables <- normalizePath(cmdArgs[3:length(cmdArgs)])

permutationTablesByGene_random <- file.path(dirname(permutationTables), "permutationsFoldChangeByGene_random.txt")
tissues <- basename(dirname(dirname(permutationTables)))

# To test
#heatmapOut <- "/scratch/users/paedugar/somaticMutationsProject/clusterMutations/cluster/cancer/driverEnrichmentPlots/individualDrivers/allDriversTissues.pdf"
#histOut <- "/scratch/users/paedugar/somaticMutationsProject/clusterMutations/cluster/cancer/driverEnrichmentPlots/individualDrivers/mean_FC_distribution.pdf"
#tissues <- list.dirs("/scratch/users/paedugar/somaticMutationsProject/clusterMutations/cluster/cancer/driverEnrichment", recursive = F)
#permutationTablesByGene <- file.path(tissues, "n6_0.0_0.5", "permutationsFoldChangeByGene.txt")
#permutationTablesByGene_random <- file.path(tissues, "n6_0.0_0.5", "permutationsFoldChangeByGene_random.txt")
#tissues <- basename(tissues)
#
#tissues <- tissues[file.exists(permutationTablesByGene)]
#permutationTables <- permutationTablesByGene[file.exists(permutationTablesByGene)]
#permutationTablesByGene_random <- permutationTablesByGene_random[file.exists(permutationTablesByGene)]

#----------------------------------#
# METHODS
#----------------------------------#

readIndividualGenesTables <- function(tissues, tables, randomTables = F) {
    
    results <- list()
    for (i in 1:length(tissues)) {
        x <- read.table(tables[i], sep = "\t", stringsAsFactors = F, header = T)
        
        if(randomTables)
            colnames(x) <- paste0("RandomGene_", 1:ncol(x))
        
        results[[i]] <- x
        results[[i]]$tissue <- tissues[i]
    }
    
    return(do.call(rbind,results))
}

plotHeatmaps <- function(signifResults, rows = "tissue", cols = "gene", vals = "realMean", heatmapsOut, zoomIn_n = 30, nCols = 20) {
    
    heatmapColor <- colorRampPalette(c("#053c93", "white", "#cec544"))(nCols)
    
    # Organizes data into 2D matrix
    dataHeat <- as.data.frame(spread(signifResults[,c(rows, cols, vals)], cols, vals))
    rownames(dataHeat) <- dataHeat$tissue
    dataHeat <- as.matrix(dataHeat[,-1])
    
    # Orders columns by mean FC and rows by mean too. NOT anymore hier clustering
    #rowOrder <- hclust(dist(dataHeat))$order
    rowOrder <- order(rowMeans(dataHeat))
    colOrder <- order(colMeans(dataHeat))
    dataHeat <- dataHeat[rowOrder,colOrder]
    
    # PLOTTING
    pdf(heatmapsOut)
    
    # Full heatmap
    maxVal <- max(abs(dataHeat))
    heatPars <- list(x = dataHeat, heatmap_gradient = heatmapColor, cexCol = 0.01, Rowv = F, Colv = F, 
                     symbreaks = T, breaks = seq(-maxVal, maxVal, length.out = nCols + 1),
                     key.title = "", key.xlab = "Mean FC (drvier/random)")
    do.call(Heatmap, heatPars)
    
    # Zoomed-in heatmaps 
    heatPars[[1]] <- dataHeat[,1:zoomIn_n]
    heatPars <- heatPars[names(heatPars) != "cexCol"]
    do.call(Heatmap, heatPars)
    
    heatPars[[1]] <- dataHeat[,(ncol(dataHeat) - zoomIn_n) : ncol(dataHeat)]
    do.call(Heatmap, heatPars)
    
    # One more heatmap with top variable genes in their FC values
    colOrder <- order(apply(dataHeat, 2, var))
    dataHeat <- dataHeat[,colOrder]
    heatPars[[1]] <- dataHeat[,(ncol(dataHeat) - zoomIn_n) : ncol(dataHeat)]
    do.call(Heatmap, heatPars)
    
    dev.off()
    
    return(invisible(NULL))
}
    
    

#----------------------------------#
# MAIN
#----------------------------------#

#----------#
# Doing plots for mean mutation fold change values of driver genes and permuted genes (same expression)

permResultsGenes <- readIndividualGenesTables(tissues, permutationTables)
permResultsGenesRandom <- readIndividualGenesTables(tissues, permutationTablesByGene_random, randomTables = T)


colnames(permResultsGenes)[-ncol(permResultsGenes)] <- ensemblToAlias(colnames(permResultsGenes)[-ncol(permResultsGenes)])

# Rearring genes
permResultsGenes <- gather(permResultsGenes, key = "gene", value = "real", -tissue)
permResultsGenesRandom <- gather(permResultsGenesRandom, key = "gene", value = "random", -tissue)
permResultsGenes$random <- permResultsGenesRandom$random

# tidy data
permResultsGenes_tidy <- gather(permResultsGenes, key = "type", value = "log2Fold_overPermutation", real, random)

# Perform tests and getting mean FCs values comparing random vs real
signifResults <- 
    permResultsGenes %>%
    group_by(tissue, gene) %>%
    summarise(pval = wilcox.test(real, random)[["p.value"]], realMean = mean(real), randomMean = mean(random), diff = realMean - randomMean) %>%
    mutate(pvalBonf = p.adjust(pval, method = "bonf"), pvalLabel = labelPvalues(pvalBonf)) %>%
    ungroup()

# Creates heatmaps of mean FC values (gene driver over randomly selected with same expression level)
plotHeatmaps(signifResults, heatmapsOut = heatmapOut)
pdf(histOut, height = 5)
hist(signifResults$realMean, col = "grey70", main = "", xlab = "mean FC (drvier/random)")
dev.off()

# Plots histograms of FCs for both random and real, ALL DATA (from all permutations)
#hist(permResultsGenes$real)
#hist(permResultsGenes$random)
#hist(signifResults$randomMean)

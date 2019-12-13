# Takes all mutation maps from tissues, it finds significant mutations based on sequencing error
# And creates individuals folders of significant and non-significant that can be used with my snakemake pipeline

library("rSubmitter")
library("purrr")
library("ggplot2")

source("../R/ggthemes.R")

# GLOBALS
BAR_WIDTH = 0.3

# PARS
projectFolder <- "/scratch/users/paedugar/somaticMutationsProject"
rootFolder <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount"

maf <- "n6_0.0_0.5"
mapFolder <- file.path(rootFolder, "map")
countFolder <- file.path(rootFolder, "count")
sequenceErrorRate <- 0.001 # qual 30
pvalueCutoff <- 0.0001

outDirClustered <- "/scratch/users/paedugar/somaticMutationsProject/significantComparison/significant"
outDirNonClustered <- "/scratch/users/paedugar/somaticMutationsProject/significantComparison/nonSignificant"

generalAnalysesPrefix <- "generalMutationAnalyses/results"
percentageVariantsInRNAsites <- file.path(outDirNonClustered, generalAnalysesPrefix, "percentageInRNA_editSites")

#-------------------------------
# METHODS

#' Given a mutation table, prints clustered mutations (distance to closest mutation < n) on file 1 and 
#' the rest of the mutations on file 2

writeSignificantMutations <- function(mutationTableFile, sequenceErrorRate, pvalue = 0.001, file_significant, file_nonSignificant) {
    
    mutationTable <- read.table(mutationTableFile, sep = "\t", header = F, stringsAsFactors = F)
    
    pvalues <- pbinom(mutationTable[,7], mutationTable[,6], sequenceErrorRate, lower.tail = F)
    
    
    # Dumping clustered mutations to one file and rest to the other
    if(sum(pvalues < pvalue) > 1)
        write.table(mutationTable[ pvalues < pvalue, ], file_significant, sep = "\t", quote = F, col.names = F, row.names = F)
    if(sum(pvalues >= pvalue) > 1)
        write.table(mutationTable[ pvalues >= pvalue, ], file_nonSignificant, sep = "\t", quote = F, col.names = F, row.names = F)

    return(NULL)
    
}

#' Print and count mutation types
#' takes fileVector, where the first element is the input map file and the second is the output count file
printCounts <- function(fileVector) {
    
    mapFile <- fileVector[1]
    outCountFile <- fileVector[2]
    
    dir.create(dirname(outCountFile), showWarnings = F, recursive = T)
    
    outMat <- matrix(0, nrow = 2, ncol = 4)
    
    flush.console()
    cat(file.info(mapFile)$size, " ", mapFile, "\n")
    
    mapTable <- tryCatch(read.delim(mapFile, sep = "\t", stringsAsFactors = F), error = function(cond) data.frame())
    
    if(nrow(mapTable) != 0) {
    
    
        frequencies <- table(paste0(mapTable[,3], mapTable[,4]))
        
        outMat[1,1] <- sum(frequencies[names(frequencies) %in% c("TA", "AT")])
        outMat[1,3] <- sum(frequencies[names(frequencies) %in% c("TG", "AC")])
        outMat[1,4] <- sum(frequencies[names(frequencies) %in% c("TC", "AG")])
        outMat[2,1] <- sum(frequencies[names(frequencies) %in% c("CA", "GT")])
        outMat[2,2] <- sum(frequencies[names(frequencies) %in% c("CT", "GA")])
        outMat[2,3] <- sum(frequencies[names(frequencies) %in% c("CG", "GC")])
        
        write.table(outMat, outCountFile, sep = "\t", row.names = F, col.names = F, quote = F)
        return(outCountFile)
    }
    
    
    return(NULL)
}
    
    
    



#-------------------------------


#-------------------------------
# MAIN



##########
# Splitting into clustered and non clustered

allTissues <- file.path(list.dirs(mapFolder, recursive = F), maf)
allFiles <- unlist(map(allTissues, function(x) file.path(x, list.files(x))))

lapply(allFiles, function(x, sequenceErrorRate, projectFolder, pvalueCutoff, outDirClustered, outDirNonClustered) {
           flush.console();cat(x, "\n")
           outCluster <- gsub(projectFolder, outDirClustered, x)
           outNonCluster <- gsub(projectFolder, outDirNonClustered, x)
           dir.create(dirname(outCluster), recursive = T,  showWarnings = F)
           dir.create(dirname(outNonCluster), recursive = T,  showWarnings = F)
           
           writeSignificantMutations(x, sequenceErrorRate, pvalueCutoff, outCluster, outNonCluster)
                   }, sequenceErrorRate = sequenceErrorRate, pvalueCutoff = pvalueCutoff, projectFolder = projectFolder, outDirClustered = outDirClustered, outDirNonClustered = outDirNonClustered)
           
           
           
    

#######
# CREATING COUNT MATRICES

allClusterFiles <- file.path(list.dirs(file.path(outDirClustered, "mutationCount", "map"), recursive = F), maf)
allClusterFiles <- c(allClusterFiles, file.path(list.dirs(file.path(outDirNonClustered, "mutationCount", "map"), recursive = F), maf))

allMapFiles <- unlist(map(allClusterFiles, ~ file.path(.x, list.files(.x)) ))
allCountFiles <- gsub("/map/", "/count/", allMapFiles)

fileList <- map2(allMapFiles, allCountFiles, ~ c(.x, .y))

allCountFiles <- superApply(fileList, printCounts, 
                       partition = "hbfraser,owners,hns,normal", mem = "2G", time = "10:00", proc = 1,
                       tasks = 200,  workingDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles",
                       packages = "")





##
# Calculate percentages of variants called in RNA edit sites
allClusterFiles <- file.path(list.dirs(file.path(outDirNonClustered, "mutationCount", "map"), recursive = F), maf)

allMapFiles <- unlist(map(allClusterFiles, ~ file.path(.x, list.files(.x)) ))
allCountFiles <- gsub("/map/", "/count/", allMapFiles)

allCountTissues <- unique(dirname(allCountFiles))

percentages <- map_dbl(allCountTissues, function(tissueFolder, originalCountFolder, maf) {
                       
                  countFiles <- file.path(tissueFolder, list.files(tissueFolder))
                  names(countFiles) <- basename(countFiles)
                  countsAll <- map_dfr(countFiles, function(countFile, originalCountFolder, maf) {
                                        
                                        sampleId <- basename(countFile)
                                        tissue <- basename(dirname(dirname(countFile)))
                                        originalCountFile <- file.path(originalCountFolder, tissue, maf, sampleId)
                                        
                                        rnaEdits <- as.matrix(read.table(countFile, sep = "\t", header = F))
                                        allVariants <- as.matrix(read.table(originalCountFile, sep = "\t", header = F))
                                        
                                        return (data.frame(rnaEdits = sum(rnaEdits), allVars = sum(allVariants)))
                                        }, originalCountFolder = originalCountFolder, maf = maf
                                        )
                  
                  return(sum(countsAll[,1]) / sum(countsAll[,2]) * 100)
                  }, originalCountFolder = countFolder, maf = maf
                  )

toPlot <- data.frame(Percent_variants_clustered = percentages, tissue = basename(dirname(allCountTissues)))
toPlot$tissue <- factor(toPlot$tissue, ordered = T, levels = toPlot$tissue[order(toPlot[,1])])

p <- ggplot(toPlot, aes(x = tissue, y = Percent_variants_clustered)) + 
geom_bar(stat = "identity") +
coord_flip() +
ylab ("Percent_variants_nonSignificant")+
ylim(c(0, 30)) +
theme_grid_x()

dir.create(percentageVariantsInRNAsites, showWarnings = F, recursive = T)
ggsave(file.path(percentageVariantsInRNAsites, "allTissues.pdf"), p, height = nrow(toPlot) * BAR_WIDTH)

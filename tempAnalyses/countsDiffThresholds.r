# Compares the mutation counts and MAFs of different Nsupport values accross different tissues
# srun -p hbfraser --time=30:00 --mem=8G --pty Rscript countsDiffThresholds.r ~/results/FraserLab/somaticMutationsProject/20180219/countsDiffThresholds/counts.pdf ~/results/FraserLab/somaticMutationsProject/20180219/countsDiffThresholds/maf.pdf

library("ggplot2")

countRootDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount"
countDir <- "count"
pileupDir <- "map"
PANEL_WIDTH <- 2
PANEL_HEIGHT <- 1.7


args <- commandArgs(T)
countFilePlot <- args[1] #"all_counts.pdf"
mafFilePlot <- args[2] #"all_maf.pdf"


#---------------------------------------------
# Methods

getTissues <- function(x) {
    
    # x in one of the count fodlers containing tissue folders
    # i.e. file.path(countRootDir, countDir)
    
    tissues <- gsub("(.+)-.+", "\\1", basename(list.dirs(x, recursive = F)))
    
    return(tissues)
}
    
getMutationPars <- function(x, tissue) {
    
    # gets a list of unique N_support values and MAFs used to 
    # to call mutations
    
    # x in one of the count fodlers containing tissue folders
    # i.e. file.path(countRootDir, countDir)
    
    features <- basename(list.dirs(file.path(x,tissue), recursive = F))
    features <- features[grep("n\\d", features)]
    
    return(features)
    
}

getCountData <- function(countRootDir, countDir, tissues, mutationFeatures) {
    
    # x the count fodler containing tissue folder and eahc of those containing counts
    # i.e. file.path(countRootDir, countDir)
    
    allCounts <- list()

    for(tissue in tissues) {
        for(mutFet in mutationFeatures){
            
            currentPath <- file.path(countRootDir, countDir, tissue, mutFet)
            mutationFiles <- file.path(countRootDir, countDir, tissue, mutFet, list.files(currentPath))
            currentCounts <- getCounts(mutationFiles)
            
            currentCounts$tissue <- tissue
            currentCounts$mutFeature <- mutFet
            
            allCounts <- c(allCounts, list(currentCounts))
        }
    }
    
    allCounts <- do.call(rbind, allCounts)
    
    return(allCounts)
}

getCounts <- function(mutationFiles) {
    
    # gets the counts as a data.frame
    
    result <- list()
    for(i in mutationFiles) {
        
        current <- read.table(i, sep = "\t", header = F, stringsAsFactors = F)
        current <- data.frame(sample = basename(i), count = sum(current), stringsAsFactors = F)
        result <- c(result, list(current))
        
    }
    
    result <- do.call(rbind, result)
    
    return(result)
    
}

plotCountDist <- function(counts, valueCol = "count", nLabel = "nInd = ") { 
    
    # Create stats matrix
    stats <- by(counts, list(counts$tissue, counts$mutFeature),
                function(x) {
                    results <- data.frame(label = paste0("mean = ", round(mean(x[,valueCol]), 2), "\n",
                                                        "sd = ", round(sd(x[,valueCol]), 2), "\n",
                                                        "max = ", round(max(x[,valueCol]), 2), "\n",
                                                        "min = ", round(min(x[,valueCol]), 2), "\n",
                                                        nLabel, nrow(x), "\n"
                                                       ),
                                          tissue =  x$tissue[1],
                                          mutFeature = x$mutFeature[1],
                                          x = Inf,
                                          y = Inf,
                                          stringsAsFactors = F
                                          )
                }
                )
    
    stats <- do.call(rbind, stats)
    
    # Plots
    counts$mutFeature <- factor(counts$mutFeature, 
                                ordered = T, 
                                levels = unique(counts$mutFeature)[order(as.numeric(gsub("n(\\d+).+", "\\1", unique(counts$mutFeature))))])
    
    stats$mutFeature <- factor(stats$mutFeature, ordered = T, levels = levels(counts$mutFeature))
    
    p <- ggplot(counts, aes_string(x = valueCol)) +
    geom_histogram() +
    facet_grid (mutFeature~tissue, scales = "free_y") +
    geom_text(aes(x = x, y = y, label = label), size = 3, data = stats, vjust = 1, hjust =1) +
    theme_bw()
    
    return(p)
}

getPileupsData <- function(countRootDir, pileupDir, tissues, mutationFeatures) {
    
    allPileups <- list()

    for(tissue in tissues) {
        for(mutFet in mutationFeatures){

            currentPath <- file.path(countRootDir, pileupDir, tissue, mutFet)
            mutationFiles <- file.path(countRootDir, pileupDir, tissue, mutFet, list.files(currentPath))
            currentPileups <- getPileups(mutationFiles)

            currentPileups$tissue <- tissue
            currentPileups$mutFeature <- mutFet

            allPileups <- c(allPileups, list(currentPileups))
        }
    }

    allPileups <- do.call(rbind, allPileups)

    return(allPileups)
}

getPileups <- function(pileupFiles) {
    
    # gets the counts as a data.frame
    
    result <- list()
    for(i in pileupFiles) {
        
        current <- read.table(i, sep = "\t", header = F, stringsAsFactors = F)
        current <- data.frame(sample = basename(i), MAF = current[,7] / current[,6], stringsAsFactors = F)
        result <- c(result, list(current))
        
    }
    
    result <- do.call(rbind, result)
    
    return(result)
    
}

    
    


#---------------------------------------------------
# Main

tissues <- getTissues(file.path(countRootDir, countDir))
mutationFeatures <- getMutationPars(file.path(countRootDir, countDir), tissues[1])

# Counts
counts <-  getCountData(countRootDir, countDir, tissues, mutationFeatures)

plotCounts <- plotCountDist(counts)
ggsave(countFilePlot, plotCounts, width = PANEL_WIDTH * length(tissues), height = PANEL_HEIGHT * length(mutationFeatures))

# MAFs
MAF <- getPileupsData(countRootDir, pileupDir, tissues, mutationFeatures)

plotMAF <- plotCountDist(MAF, "MAF", "nMut = ")
ggsave(mafFilePlot, plotMAF, width = PANEL_WIDTH * length(tissues), height = PANEL_HEIGHT * length(mutationFeatures))


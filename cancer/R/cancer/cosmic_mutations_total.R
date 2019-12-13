# Calculates how many of the ALL the cosmic mutations were observed at least once in this study
library("ggplot2")
source("../../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    mutationDir <- cmdArgs[1]
    maf <- cmdArgs[2]
    
    
    mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/map"
    maf <- "n6_0.0_0.7"
    
    counts <- read_cosmic_counts(mutationDir, maf)
    
    
}

read_cosmic_counts <- function (mutationDir, maf) {
    
    results <- list()
    i <- 1 
    for(tissue in list.dirs(mutationDir, recursive = F)) {
        if (grepl("EXO", tissue))
            next
        for(current_sample in list.files(file.path(tissue, maf))) {
            current <- read.table(file.path(tissue, maf, current_sample), sep = "\t", stringsAsFactors = F, header = F)
            current <- current[current[,5] == current[,10], 1:5]
            #if(!exists(results)) {
            #    results <- current
            #} else {
            #    results <- rbind(results, current)
            #    results <- unique(results)
            #}
            results[[i]] <- current
            i <- i + 1
        }
    }
    
    results <- do.call(rbind, results)
    #results <- unique(results)
    
    
    return(results)
    
}

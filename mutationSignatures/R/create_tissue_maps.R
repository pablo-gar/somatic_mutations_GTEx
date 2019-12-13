main <- function(cmdArgs = commandArgs(T)) {
    
    out_map_file <- cmdArgs[1]
    mutation_map_files <- cmdArgs[-1]
    
    mutation_map <- readMutationMaps(mutation_map_files)
    
    write.table(mutation_map, out_map_file, sep = "\t", quote = F, col.names = F, row.names = F)
    
}

readMutationMaps <- function(mutationFiles) {
    
    mutationsAll <- list()
    for(i in 1:length(mutationFiles)) {
        
        mutationFile <- mutationFiles[i]
        
        currentSample <- gsub(".txt", "", basename(mutationFile))
        
        mutations <- read.table(mutationFile, sep = "\t", stringsAsFactors = F)
        
        colnames(mutations) <- c("chr", "pos", "ref", "mut", "context", "coverage", "nSupport")
        mutations$sample <- currentSample
        mutationsAll[[i]] <- mutations
        closeAllConnections()

    }

    mutationsAll <- do.call(rbind, mutationsAll)
    
    mutationsAll <- mutationsAll[,c("sample", "chr", "pos", "ref", "mut")]
    
    return(mutationsAll)
}

main()

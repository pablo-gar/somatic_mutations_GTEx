library("ggplot2")
library("rSubmitter")

sourceForCluster <- "~/scripts/FraserLab/somaticMutationsProject/mutationSignatures/R/mutationSignatures/source/NMF.class.R"

readMutations <- function(x, contextLength, ignoreN) {
    
    # Reads a file with columns as Mutation Type, Context, Ind 1, ... , Ind n
    # Merges rows with identical contexts depending on the context length


    #Getting number of columns
    con <- file(x, open = "r")
    columnsFile <- length(unlist(strsplit(readLines(con, n = 1), "\t")))
    close(con)


    x <- read.delim(x, header = T, sep = "\t", stringsAsFactors = F, colClasses = c("character", "character", rep("numeric", columnsFile - 2)) )

    if (ignoreN)
        x <- x[!grepl("N", x[,2]), ]

    currentLength <-  nchar(x[1,2])

    if (contextLength > currentLength)
        stop("context length specified greater than currently availabe")

    if (contextLength < currentLength) {
        lengthDiff <- (currentLength - contextLength) / 2
        x[,2] <- substr(x[,2], lengthDiff + 1, currentLength - lengthDiff)

        # Put same context together
        x <- do.call(rbind, by(x, list(x[,1], x[,2]), function(x) {
                                   x[1,3:ncol(x)] <- colSums(x[,3:ncol(x)])
                                   return (x[1,])
                                    }
                                )
                    )

    }


    # Reordering and renaming
    x <- x[order(x[,1], x[,2]),]
    rownames(x) <- 1:nrow(x)
    return(x)
}

readContextCounts <- function(mutationMat, oligoCountsFile) {

    # Reads the oligo count file into a vector, makes sure that oligo length
    # is equal to context length in mutation file.
    # Reorders oligo count file such that each position correspond to
    # each row in the mutation matrix

    oligoCounts <- read.table(oligoCountsFile, sep = "\t", stringsAsFactors = F, header = F)
    rownames(oligoCounts) <- oligoCounts[,1]
    

    if (nchar(oligoCounts[1,1]) != nchar(mutationMat[1,2]))
        stop ("context length is not the same as oligonucleotide length")

    if (sum(!mutationMat[,2] %in% rownames(oligoCounts)) > 0)
        stop("some context types are missing in oligo count table")

    result <- oligoCounts[mutationMat[,2], 2]
    names(result) <- paste0(mutationMat[,1], ".", mutationMat[,2])

    return(result)

}

estimateNsign <- function (mutationMat, normalizingVector, tryNsignatures, nRunsPerNMF, totalSimulations, useFreq, nCores, clusterWorkingDir) {
    
    similarity <- rep(0, length(tryNsignatures))
    fError <- rep(0, length(tryNsignatures))
    
    # Trying different number of signatures
    for (i  in 1:length(tryNsignatures)) {
        nSign <- tryNsignatures[i]
        
        flush.console()
        cat("Working with ", nSign, " signatures\n")
        
        
        # Running simulations
        allSimulations <- "failed"
        tries <- 1
        while(is.character(allSimulations)) {
        cat("Running", tries, " try\n")
        allSimulations <- tryCatch(superApply(1:totalSimulations, workingDir = clusterWorkingDir, tasks = round(totalSimulations/2), mem = "16G", time = "06:00:00",
                                        sources = sourceForCluster, packages = "", 
                                        partition = "hbfraser,normal,hns",
                                        mutationMat = mutationMat, normalizingVector = normalizingVector,
                                        useFreq = useFreq, reducedPercentage = 0.01, nRuns = nRunsPerNMF,
                                        nSign = nSign, nCores = nCores,
                                        FUN = function(X, mutationMat, normalizingVector, useFreq, reducedPercentage, nRuns, nSign, nCores) {
                                            result <- NMFmontecarloList$new(mutationMat = mutationMat, normalizingVector = normalizingVector, useFreq = useFreq, 
                                                                  reducedPercentage = 0.01, nRuns = nRuns, nSign = nSign, nCores = nCores)
                                            return(result)
                                        }),
                                    error = function(e) { print(e); "failed" }
                                    )
        
        if (tries >= 4)
            stop("Reached 4 tries stopping now")
        tries <- tries + 1
        }
        
        
        
        
        # Merging simulations
        nmfList <- allSimulations[[1]]$clone()
        for(j in 2:(length(allSimulations))) {
            nmfList$c(allSimulations[[j]])
        }
        
        # Calcualte and store similarities of a given signature number
        similarity[i] <- nmfList$calculateSimilarity()
        fError[i] <- nmfList$calculateFrobError()
    }
    
    # Select the number of signatures where similarity and error difference is the highest
    maxDiff <- rank(-similarity) + rank(fError)
    maxDiff <- min(maxDiff) == maxDiff
    bestNsign <- tryNsignatures[maxDiff][1]
    
    # Plot results
    toPlot <- data.frame(nSignatures = rep(tryNsignatures, 2), value = c(similarity, fError), 
                         metric = rep(c("similarity", "fError"), each = length(tryNsignatures)), 
                         selected = rep(maxDiff, 2))
    
    p <- ggplot(toPlot, aes(x = nSignatures, y = value))+
    geom_point(aes(colour = selected)) +
    facet_grid(metric~., scale = "free_y") +
    theme_bw()
                
    
    return(list(nSign = bestNsign, plot = p))
    
    
}



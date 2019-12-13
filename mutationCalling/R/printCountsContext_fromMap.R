#' Takes a mutation map matrix and prints two other useful matrices, a count matrix and a 
#' context matrix.
#'
#' @param context length, hast to be an odd number
#' @param file path to map matrix: tab-separated file with columns chr, coordinate, ref base, alt base, sequence context, coverage, alt counts
#' @param file path to output count matrix, first row is ref base T, second row ref is C, columns are alt alleles, A, T, G, C,
#' @param file path to output context matrix with columns, mutation type, context, count
#' 

#' Print and count mutation types
#' takes fileVector, where the first element is the input map file and the second is the output count file

library("gtools")

main <- function(cmdArgs = commandArgs(T)) {
    
    cntxtLength <- as.integer(cmdArgs[1])
    mapFile <- cmdArgs[2]
    outCountFile <- cmdArgs[3]
    outCntxtFile <- cmdArgs[4]
    
    #Debugging
    #cntxtLength <- 5
    #mapFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCountIncludeMajor/map/Adipose_Subcutaneous/n6_0_1/SRR1069097.txt"
    #mapFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Whole_Blood_EXO/n6_0.0_0.7/SRR2165971.txt"
    #mapFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Spleen/n6_0.0_0.7/SRR1102296.txt"
    #outCountFile <- "outCount.txt"
    #outCntxtFile <- "outCntxt.txt"
    
    if(!file.exists(mapFile))
        stop("File does not exist:", mapFile)
    
    mapTable <- readMapTable(mapFile, cntxtLength)
    
    # Correct bases where the reference flipped or is incorrect at all
    mapTable <- correct_bases(mapTable, cntxtLength)
    write.table(mapTable, mapFile, sep = "\t", col.names = F, row.names = F, quote = F)
    
    printCount(mapTable, outCountFile)
    printCntxt(cntxtLength, mapTable, outCntxtFile)
    
}

correct_bases <- function(mapTable, cntxtLength) {
    middleBase <- ceiling(cntxtLength/2)
    refBases <- substring(mapTable$V5, middleBase, middleBase)
    # Remove the complete incorrect ones
    mapTable <- mapTable[!(mapTable$V3 !=  refBases & mapTable$V4 !=  refBases),]
    
    # Correct others
    refBases <- substring(mapTable$V5, middleBase, middleBase)
    toCorrect <- mapTable$V4 == refBases
    
    if(any(toCorrect)) {
        V4 <- mapTable[toCorrect, "V4"]
        mapTable[toCorrect, "V4"] <- mapTable[toCorrect, "V3"]
        mapTable[toCorrect, "V3"] <- V4
    }
    
    # Eliminate non SNVS
    mapTable <- mapTable[nchar(mapTable$V3) == 1 & nchar(mapTable$V4) == 1,]
    return(mapTable)
    
}
    
readMapTable <- function(mapFile, cntxtLength) {
    
    mapTable <- tryCatch(read.delim(mapFile, sep = "\t", header = F, stringsAsFactors = F, 
                                    colClasses = c("character", "numeric", "character", "character", "character", "numeric", "numeric" ))
                                    , error = function(err) data.frame())
    
    if(nrow(mapTable) == 0)
        return(mapTable)
    
    if(nchar(cntxtLength) > nchar(mapTable[1,5]))
        stop("Specified context length is greater than context length in map file")
    
    # If desired context length is less than specified, shrink the context
    if(cntxtLength < nchar(mapTable[1,5])) {
        diffL <- (nchar(mapTable[1,5]) - cntxtLength) / 2
        mapTable[,5] <- substr(mapTable[,5], 1 + diffL,  nchar(mapTable[1,5]) - diffL)
    }
    
    return(mapTable)
    
}

printCount <- function(mapTable, outCountFile) {
    
    
    outMat <- matrix(0, nrow = 2, ncol = 4)
    
    if(nrow(mapTable) != 0) {
    
        frequencies <- table(paste0(mapTable[,3], mapTable[,4]))
        
        outMat[1,1] <- sum(frequencies[names(frequencies) %in% c("TA", "AT")])
        outMat[1,3] <- sum(frequencies[names(frequencies) %in% c("TG", "AC")])
        outMat[1,4] <- sum(frequencies[names(frequencies) %in% c("TC", "AG")])
        outMat[2,1] <- sum(frequencies[names(frequencies) %in% c("CA", "GT")])
        outMat[2,2] <- sum(frequencies[names(frequencies) %in% c("CT", "GA")])
        outMat[2,3] <- sum(frequencies[names(frequencies) %in% c("CG", "GC")])
        
    }
    
    write.table(outMat, outCountFile, sep = "\t", row.names = F, col.names = F, quote = F)
    
}

printCntxt <- function(cntxtLength, mapTable, outCntxtFile) {
    
    middleBase <- ceiling(cntxtLength/2)
    translate <- setNames(c("A", "T", "G", "C", "N"), c("T", "A", "C", "G", "N"))
    from <- setNames(c("0", "1"), c("T", "C"))
    to <- setNames(c("0", "1", "2", "3"), c("A", "T", "G", "C"))
    
    cntxtMat <- createCntxtMat(cntxtLength)
    
    # Get ready mutations
    for(i in seq_len(nrow(mapTable))) {
        ref = mapTable[i,3]
        alt = mapTable[i,4]
        cntxt = mapTable[i,5]
        
        
        if(ref %in% c("A", "G")) {
            ref <- translate[ref]
            alt <- translate[alt]
            cntxt <- paste0(substring(cntxt, 1, middleBase - 1), ref, substring(cntxt, middleBase + 1))
        } 
        
        mutType <- paste0(from[ref], to[alt])
        codeMut <- paste0(mutType, cntxt)
        cntxtMat[codeMut, "count"] <- cntxtMat[codeMut, "count"] + 1
    }
    
    if(any(is.na(cntxtMat$cntxt))) {
        print(cntxtMat[is.na(cntxtMat$cntxt),])
        stop("Something went wrong, a mutation type was not expected")
    }
    
    write.table(cntxtMat, outCntxtFile, sep = "\t", col.names = F, row.names = F, quote = F)
    
}

createCntxtMat <- function(cntxtLength) {
    
    middleBase <- ceiling(cntxtLength/2)
    cntxtAll <- permutations(n = 5, r = cntxtLength, v = c("A", "T", "G", "C", "N"), repeats.allowed = T)
    cntxtAll <- cntxtAll[cntxtAll[ ,middleBase] != "N", ]
    cntxtAll[cntxtAll[ ,middleBase] == "A" ,middleBase] <- "T"
    cntxtAll[cntxtAll[ ,middleBase] == "G" ,middleBase] <- "C"
    
    cntxtT <- cntxtAll[cntxtAll[ ,middleBase] == "T", ]
    cntxtC <- cntxtAll[cntxtAll[ ,middleBase] == "C", ]
    
    cntxtT <- apply(cntxtT, 1, paste, collapse = "")
    cntxtC <- apply(cntxtC, 1, paste, collapse = "")
    
    cntxtT <- unique(cntxtT)
    cntxtC <- unique(cntxtC)
    
    cntxtMat <- data.frame(mutationType = rep(c("00", "02", "03", "10", "11", "12"), each = length(cntxtT)), 
                           cntxt = c(rep(cntxtT, 3), rep(cntxtC, 3)), 
                           count = 0)
    rownames(cntxtMat) <- paste0(cntxtMat$mutationType, cntxtMat$cntxt)
    
    return(cntxtMat)
}


main()

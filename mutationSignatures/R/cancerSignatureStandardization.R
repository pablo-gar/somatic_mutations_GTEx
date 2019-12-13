## Takes the mutation signature file from cosmic and returns 
## the standard W matrix as it is return from NMF



args <- commandArgs(T)
cancerSignaturesFile <- args[1]
outFile <- args[2]

readCancerSignatures <- function(x) {
    
    x <- read.delim(x, sep = "\t", stringsAsFactors = F, header = T)
    rownames(x) <- paste0(toNumericMutation(x[,1]), ".", x[,2])
    x <- x[,-(1:3)]
    
    return(x[,1:30])
}

#' Converts mutation from string code to numeric, e.g. T>A = 00
toNumericMutation <- function(x) {
    
    MUTATION_KEY <- c("T>A" = "00", "T>T" = "01", "T>G" = "02", "T>C" = "03",
                      "C>A" = "10", "C>T" = "11", "C>G" = "12", "C>C" = "13")
    
    return(MUTATION_KEY[x])
    
}

cancerSignatures <- readCancerSignatures(cancerSignaturesFile)
write.table(cancerSignatures, outFile, sep = "\t", col.names = F, row.names = T, quote = F)

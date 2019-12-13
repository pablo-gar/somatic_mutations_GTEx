# Calculates shannon diversity for all samples in gtex
# Usage
# Rscript get_transcriptome_diversity.R out_table.txt

source("../../R/gtex.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    out_file <- cmdArgs[1]
    
    out_table <- get_shannon()
    write.table(out_table, out_file, sep = "\t", quote = F, row.names = F, col.names = T)
}


# calculates the shanon entropy of a vector
shannonEntropy <- function(x, percentage = T) {
    
    x <- x[ x != 0 & !is.na(x) ]
    
    shan <- -sum( (x/sum(x)) * log2(x/sum(x)) )
    
    if(percentage)
        shan <- shan / log2(length(x))
    
    return(shan)
}

# Calcautes the shanon entropy across columns in a matrix in an efficent way
# IT DOES NOT ELIMINATE 0s
shannonEntropyMatrix <- function(x) {
    
    small <- 0.000001
    
    # Doing frequencies by column
    x <- x / matrix(colSums(x), nrow = nrow(x), ncol = ncol(x), byrow=T)
    
    shan <- -colSums((x + small) * log2(x + small))
    
    return(shan)
}

get_shannon <- function(ignoreZero = T) {

    geneExp <- readAllGeneReadCounts()
    
    if(ncol(geneExp) < 2 | nrow(geneExp) < 2)
        return(NULL)
    
    geneExp <- geneExp[,-1]
    
    if(ignoreZero) {
        shan <- apply(as.matrix(geneExp), 2, shannonEntropy)
    } else {
        shan <- shannonEntropyMatrix(as.matrix(geneExp))
    }
    

    shan <- data.frame(gtexId = names(shan), transcriptome_diversity = shan)
    
    return(shan)
    
}

main()

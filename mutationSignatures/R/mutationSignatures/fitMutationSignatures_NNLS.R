#
# Given a mutation count matrix with rows being different mutation types
# and columns different individuals, a cancer signature W matrix (cosmic) and 
# a presence matrix (signatures in cancer types) it performs nnls to find
# the weight of each signature in a tumor sample
#
# It can normalize by a context count table and calculate frequencies
#
# Author: Pablo E. Garcia-Nieto
# 
# Date: 6/1/2018
# 
# Params
#   @args - mutation context length
#   @args - should I normalize by context absolute counts? (TRUE/FALSE)
#   @args - working dir to store temporary cluster files
#   @args - context count matrix
#   @args - cancer signature file from cosmic
#   @args - boolean matrix with columns as cancer types and rows as signatures indicating the presence of a signature in a cancer type
#   @args - mutation count matrix
#   @args - output file name (H matrix) 

# Rscript fitMutationSignatures_NNLS.R 2 F /scratch/users/paedugar/tcga_mutationSignatureProject/clusterFiles/ /scratch/users/paedugar/tcga_mutationSignatureProject/auxiliaryFiles/GTF_gtex_oligoCounts_3.txt /scratch/users/paedugar/tcga_mutationSignatureProject/auxiliaryFiles/cancerSignatures.txt  /scratch/users/paedugar/tcga_mutationSignatureProject/auxiliaryFiles/cancer_signature_prescence_TCGA.txt /scratch/users/paedugar/tcga_mutationSignatureProject/mutationSignatures/mutationMatrix/BLCA.txt /scratch/users/paedugar/tcga_mutationSignatureProject/mutationSignatures/mutationSignaturesNMF/BLCA_new.txt

library("nnls")
library("purrr")
library("gplots")

##-------------------------------------------
## Args and globals

args <- commandArgs(T)
contextLength <- as.integer(args[1])
normalizeContext <- as.logical(args[2])
oligoCountsFile <- args[3]
cancerSignaturesFile <- args[4]
cancerSignaturesPresenceFile <- args[5]
mutationFile <- args[6]
outFileH <- args[7]

ignoreN <- TRUE

# Debugging
#contextLength <- 3
#normalizeContext <- T
#mutationFile <-  "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationMatrix/Colon_Transverse-n6_0.0_0.7.txt"
#oligoCountsFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n3"
#cancerSignaturesFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/consensusMutationSignatures/n6_0.0_0.7_consensusSignatures.txt"
#outFileH <- "~/H.txt"
#cancerSignaturesPresenceFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/consensusMutationSignatures/n6_0.0_0.7_presenceSignatures.txt"




##-------------------------------------------
## METHODS

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

readCancerSignatures <- function(x) {
    
    x <- read.delim(x, sep = "\t", stringsAsFactors = F, header = F, row.names = 1)
    return(x)
}


##-------------------------------------------


    
##-------------------------------------------
## MAIN

#####
# Reading files
mutationType <-  basename(mutationFile)

mutationMat <- readMutations(mutationFile, contextLength, ignoreN)
cancerSignatures <- readCancerSignatures(cancerSignaturesFile)

cancerSignaturesPresence <- read.table(cancerSignaturesPresenceFile, sep = "\t", stringsAsFactors = F, header = T, check.names = F)
cancerSignaturesPresence <- cancerSignaturesPresence[,mutationType]

if(normalizeContext) {
    contextCounts <- readContextCounts(mutationMat, oligoCountsFile)
} else {
    contextCounts <- NULL
}


# Create a numeric only mutation matrix
rownames(mutationMat) <- paste0(mutationMat[,1], ".", mutationMat[,2])
mutationMat <- mutationMat[,-(1:2)]

# throwing away individuals with no info
mutationMat <- mutationMat[,colSums(mutationMat) != 0]

#cancerSignatures <- as.matrix(cancerSignatures[rownames(mutationMat), cancerSignaturesPresence, drop = F])
cancerSignatures <- as.matrix(cancerSignatures[rownames(mutationMat),, drop = F])

if(normalizeContext) 
    mutationMat <- mutationMat / contextCounts

## DONE
####

#####
## Running NNLS

H <- map_dfr(as.list(mutationMat), function(x, cancerSignatures){
                       coefs <- coef(nnls(cancerSignatures, x))
                       normCoefs <- coefs / sum(coefs)
                       return(normCoefs)
                      }, cancerSignatures = cancerSignatures)
H$signature <- colnames(cancerSignatures)
H <- H[, c(ncol(H), 2:(ncol(H)-1))]

# Eliminating the "Signature" text
H$signature <- as.numeric(gsub("V", "", H$signature)) - 1

toPlot <- as.data.frame(H)
rownames(toPlot) <- toPlot[,1]
toPlot <- toPlot[,-1]

pdf(paste0(outFileH, ".heatmap.pdf"))
heatmap.2(as.matrix(toPlot), col = colorRampPalette(c("White", "Blue"))(10), Colv = F, Rowv = F, trace = "none")
dev.off()

# Outputting files
write.table(H, outFileH, sep = "\t", quote = F, row.names = F)

## Given a mutation signature table (table W from NMF) 
## List the highest similar cancer signature to all the signatures found in the table
##
## Similarity is based on cosine similarity, and a p-value is calculated based on a permutation
## strategy, whereby rows of the NMF table are permuted n times and cosine similarity is recalculated
## 
## Author: Pablo E. Garcia-Nieto
## Date: 2/28/2018
##
## Usage:
##  Rscript signatureCancerSimilarity.R cancerSignatures.txt Skin_Sun_Exposed_Lower_leg-n4_0.0_0.5_W.txt 100000
##
## @return Prints to screen a tab-separated table with the following columns
##    1. name of W file
##    2. original_signature
##    3. cancer signature
##    4. cosine similarity
##    5. permutation number
##    5. cosine similarity
##    6. pvalue



args <- commandArgs(T)
cancerSignaturesFile <- args[1]
Wfile <- args[2]
permutations <- as.numeric(args[3])

#--------------------------------
# METHODS

readCancerSignatures <- function(x) {
    
    x <- read.delim(x, sep = "\t", stringsAsFactors = F, header = F, row.names = 1)
    return(x)
    
}

#' Converts mutation from string code to numeric, e.g. T>A = 00
toNumericMutation <- function(x) {
    
    MUTATION_KEY <- c("T>A" = "00", "T>T" = "01", "T>G" = "02", "T>C" = "03",
                      "C>A" = "10", "C>T" = "11", "C>G" = "12", "C>C" = "13")
    
    return(MUTATION_KEY[x])
    
}

#' Calculates cosine similarities between the corresponding columns
#' of two matrices
cosineSimilarity <- function(a, b) {
    
    stopifnot(dim(a) == dim(b))
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b))
    
    colSums(a * b) / (sqrt(colSums(a^2)) * sqrt(colSums(b^2)))
}

#' Performs cosine similarity and permutes a
cosineSimilarityPermute <- function(a, b) {
    
    a <- a[sample(nrow(a), nrow(a)),]
    return(cosineSimilarity(a,b))
    
}

cosineSimilarityPermuteLapply <- function(x,a,b){
    return(cosineSimilarityPermute(a,b))
}
    



#--------------------------------




#--------------------------------
# MAIN

W <- read.table(Wfile, sep = "\t", stringsAsFactors = F, header = F, row.names = 1)
cancerSignatures <- readCancerSignatures(cancerSignaturesFile)
cancerSignatures <- cancerSignatures[rownames(W), ]

# Get most similar signature to all signatures in W
mostSimilar <- rep(-1, ncol(W))
names(mostSimilar) <- "none"
for(j in 1:ncol(cancerSignatures)) {
    
    sign <- cancerSignatures[,j]
    signName <- colnames(cancerSignatures)[j]
    
    cancerMatrix <- matrix(rep(sign, ncol(W)), ncol = ncol(W), byrow = F)
    cosineResult <- cosineSimilarity(as.matrix(W), cancerMatrix)
    
    # Assesing whether this signature is more similar to previous ones
    isBetter <- mostSimilar < cosineResult
    if(any(isBetter)) {
        mostSimilar[isBetter] <- cosineResult[cosineResult >= mostSimilar]
        names(mostSimilar)[isBetter] <- signName
    }
    
}

# Permute and see calculate p values based on how many time I observed a better match in the search
permuted_i <- lapply(logical(permutations), function(x,y) return(sample(nrow(y), replace = F)), y = W)
pvalue <- rep(1, length(mostSimilar))
for(i in seq_along(mostSimilar)) {
    permuted_W <- as.matrix(do.call(cbind, lapply(permuted_i, function(x,y) y[x], y = W[,i])))
    better <- 0
    total <- 0
    for(j in 1:ncol(cancerSignatures)) {
        
        sign <- cancerSignatures[,j]
        signName <- colnames(cancerSignatures)[j]
        
        cancerMatrix <- matrix(rep(sign, ncol(permuted_W)), ncol = ncol(permuted_W), byrow = F)
        cosineResult_perm <- cosineSimilarity(permuted_W, cancerMatrix)
        
        # Assesing how many time is better
        better <- better + sum(cosineResult_perm >= mostSimilar[i])
        total <- total + length(cosineResult_perm)
    }
    
    pvalue[i] <- better/total
    
}




# Finally prints to screen a table with the following columns
cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "file", "signature", "cancerSignature", "permutations", "pvalue", "total_perm", "cosine"))
for(j in 1:ncol(W)) 
    cat(sprintf("%s\t%.0f\t%s\t%.0f\t%f\t%.0f\t%f\n", basename(Wfile), j, names(mostSimilar)[j], permutations, pvalue[j], total, mostSimilar[j]))

#--------------------------------

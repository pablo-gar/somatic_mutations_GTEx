#
# Given a mutation count matrix with rows being different mutation types
# and columns different individuals it performs Non-negative Matrix Factorization to find
# r number of mutations signatures
#
# It can normalize by a context count table and calculate frequencies
#
# Applies Burnet method and uses a unifomr random seeding.
# Outputs two matrices, W and H. Where
#
# X = WH
# 
# W is the mutation signature matrix(mutation types x mutation signatures)
# H is the contriubtion matrix(individuals x mutation signatures)
#
#
# Author: Pablo E. Garcia-Nieto
# 
# Date: 11/27/2017
# 
# Params
#   @args - number of independent runs
#   @args - number of signatures (ignore if discover signatures set to TRUE)
#   @args - should I discover the number of signatures by myself? (TRUE/FALSE)
#   @args - mutation context length
#   @args - should I use frequencies over individuals? (TRUE/FALSE)
#   @args - should I normilize by context absolute counts? (TRUE/FALSE)
#   @args - number of cores to use
#   @args - working dir to store temporary cluster files
#   @args - mutation count matrix
#   @args - context count matrix
#   @args - output file name (W matrix) 
#   @args - output file name (H matrix) 
#   @args - output file name (plot error and reproducibility) 

# Rscript mutationSignaturesNMF.r 10000 5 T 3 T T 16 /scratch/users/paedugar/deletemeSuperApply/ /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/mutationMatrixTest.txt /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n3 ~/W.txt ~/H.txt ~/plot.pdf
library("NMF")
library("R6")

##-------------------------------------------
## Args and globals

args <- commandArgs(T)
nRuns <- as.integer(args[1])
nSign <- as.integer(args[2])
discoverNsign <- as.logical(args[3])
contextLength <- as.integer(args[4])
useFreq <- as.logical(args[5])
normalizeContext <- as.logical(args[6])
nCores <- as.integer(args[7])
clusterWorkingDir <- args[8]
oligoCountsFile <- args[9]
mutationFile <- args[10]
outFileW <- args[11]
outFileH <- args[12]
outFilePlot <- args[13]

ignoreN <- TRUE
nRunsForNSignatures <- 250 # Number of runs within each NMF for N signatures estimation
totalSimulations <- 500 # Number of simulated mutation matrices for N signatures estimation
signaturesN <- 2:6

#ignoreN <- TRUE
#nRunsForNSignatures <- 200 # Number of runs within each NMF for N signatures estimation
#totalSimulations <- 100 # Number of simulated mutation matrices for N signatures estimation
#signaturesN <- 2:6

# Debugging
#nSignatures <- 4
#nRuns <- 5000
#nSign <- 4
#discoverNsign <- F
#contextLength <- 3
#useFreq <- F
#normalizeContext <- T
#nCores <- 16
#mutationFile <-  "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationMatrix/Lung-n6_0.0_0.7.txt"
#oligoCountsFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n3"
#ignoreN <- TRUE
#nRunsForNSignatures <- 200 # Number of runs within each NMF for N signatures estimation
#totalSimulations <- 100 # Number of simulated mutation matrices for N signatures estimation
#signaturesN <- 2:3
#clusterWorkingDir <- "/scratch/users/paedugar/deletemeSuperApply/"
#outFileW <- "~/W_Lung.txt"
#outFileH <- "~/H_Lung.txt"
#outFilePlot <- "~/Lung_plot.pdf"



sourceDir <- function (path, pattern = "\\.[rR]$", env = NULL, chdir = TRUE) {
    files <- sort(dir(path, pattern, full.names = TRUE))
    lapply(files, source, chdir = chdir)
    return(NULL)
}

    
##-------------------------------------------
## MAIN
sourceDir("~/scripts/FraserLab/somaticMutationsProject/mutationSignatures/R/mutationSignatures/source/")

#####
# Reading files
mutationMat <- readMutations(mutationFile, contextLength, ignoreN)
if(normalizeContext) {
    contextCounts <- readContextCounts(mutationMat, oligoCountsFile)
} else {
    contextCounts <- NULL
}

# Create a numeric only mutation matrix
rownames(mutationMat) <- paste0(mutationMat[,1], ".", mutationMat[,2])
mutationMat <- mutationMat[,-(1:2)]

# throwing away individuals with no info and mutation types with no info
mutationMat <- mutationMat[,colSums(mutationMat) != 0]
mutationTypesMissing <- rownames(mutationMat)[rowSums(mutationMat) == 0]

## DONE
####

#####
## Running NMF

if (discoverNsign) {
    nSign <- estimateNsign (mutationMat = mutationMat, normalizingVector = contextCounts, tryNsignatures = signaturesN, 
                               nRunsPerNMF = nRunsForNSignatures, totalSimulations = totalSimulations,
                               useFreq = useFreq, nCores = 1,  clusterWorkingDir = clusterWorkingDir)
    ggsave(outFilePlot, nSign$plot)
} else {
    nSign <- list(nSign = nSign)
}
    
# Aplying NMF
nmfResult <- NMF$new(mutationMat, normilizingVector = contextCounts, useFreq = useFreq, nRuns = nRuns, nSign = nSign$nSign, nCores = nCores)
nmfResult$nmf()
W <- nmfResult$W
H <- nmfResult$H


# Outputting files
write.table(W, outFileW, sep = "\t", quote = F, col.names = F)
write.table(H, outFileH, sep = "\t", quote = F, row.names = F)

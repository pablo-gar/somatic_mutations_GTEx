# Creates table of correlations between mutation signature value and age of 
# individuals
#
# Takes H matrices (weight matrices) from NMF mutation signature algorithm
#
# Author: Pablo E. Garcia-Nieto

# Usage
# Rscript ageCorrelation.R ageCorrelationHistogram.pdf ageMutationScatterPlot.pdf Hmat1.txt [Hmat2.txt ...]
# Rscript ageCorrelation.R FALSE /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/ageCorrelation/ageCorrelationHistogram_n6.pdf /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/ageCorrelation/ageMutationScatterPlot_n6.pdf /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/*n6*H.txt > /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/ageCorrelation/ageCorrelationTable_n6.txt
# Rscript ageCorrelation.R TRUE /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/ageCorrelation_consensus/ageCorrelationHistogram_n6.pdf /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/ageCorrelation_consensus/ageMutationScatterPlot_n6.pdf /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/*n6*H_consensus.txt > /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/ageCorrelation_consensus/ageCorrelationTable_n6.txt

# @return
#   1. Prints a correlation table to screen, ordered in increasing order of correlation coefficient
#   2. Saves 2 plots, a histogram of corrlation coefficient and scatter plots for correlations greater than 0.3

library(ggplot2)
source("~/scripts/FraserLab/somaticMutationsProject/R/gtex.R", chdir = T)
source("~/scripts/FraserLab/somaticMutationsProject/R/ggthemes.R", chdir = T)
source("~/scripts/FraserLab/somaticMutationsProject/R/plots.R", chdir = T)
       
# Global
PANEL_WIDTH = 4
PANEL_HEIGHT = 4
#SCATTER_GREATER_THAN = 0.25
#SCATTER_SMALLER_THAN = -SCATTER_GREATER_THAN
BONFERRONI_SIGNIF = 0.05

# Gathering arguments
args <- commandArgs(T)
signatureColumn <- as.logical(args[1])
outHistFile <- args[2]
outScatterFile <- args[3]
args <- args[-(1:3)]

dir.create(basename(outHistFile), recursive = T)


# Performing correlations
fileName <- vector()
signature <- vector()
corResults <- vector()
pvalue <- vector()
                    
for(i in 1:length(args)) {
    
    Hfile <- args[i]
    
    # Reading table
    Hmat <- read.table(Hfile, header = T, sep = "\t", stringsAsFactors = T)
    if(signatureColumn) {
        signatureNames <- as.character(Hmat[,1])
        Hmat <- Hmat[,-1]
    } else {
        signatureNames <- as.character(1:nrow(Hmat))
    }
    rownames(Hmat) <- signatureNames
    
        
    age <- readMetadataGTEX("AGE")[sraToGtex(colnames(Hmat)), , drop = F ]
    colnames(Hmat) <- sraToGtex(colnames(Hmat))
    age <- age[rownames(age) %in% colnames(Hmat),,drop=F]
    
    # Calculate frequencies
    #Hmat <- t(t(Hmat) / (colSums(Hmat)))
    
    fileName <- c(fileName, rep(Hfile, nrow(Hmat)))
    flush.console()
    for(r in 1:nrow(Hmat)) { 
        currentCorResult <- cor.test(as.numeric(Hmat[r,rownames(age)]), age[,1])
        
        signature <-  c(signature, signatureNames[r])
        corResults <- c(corResults, as.numeric(currentCorResult$estimate))
        pvalue <- c(pvalue, currentCorResult$p.value)
    }
    

}

results <- data.frame(file = fileName, signature = signature, cor = corResults, pvalue = pvalue, fdr = 0, stringsAsFactors = F)
results$fdr <- p.adjust(results$pvalue)
results <- results[order(results$cor),]
# Print results to screen
cat(paste0(colnames(results), collapse = "\t"), "\n")
for(r in 1:nrow(results))
    cat(paste0(results[r,], collapse = "\t"), "\n")


# Creating histogram of correlations
pHist <- ggplot(results, aes(x = signif(cor, 2))) +
         geom_histogram() +
         xlab("Person coefficient") +
         theme_bw()
    
# Creating scatter plots for correlations >0.3
#signifCor <- results[results$cor >= SCATTER_GREATER_THAN | results$cor <= SCATTER_SMALLER_THAN, ]
signifCor <- results[results$fdr < BONFERRONI_SIGNIF, ]

age_all <- vector()
signature <- vector()
value <- vector()

for(i in 1:nrow(signifCor)) {
    
    Hfile <- signifCor[i,"file"]
    
    # read Hmats
    Hmat <- read.table(Hfile, header = T, sep = "\t", stringsAsFactors = T)
    if(signatureColumn) {
        signatureNames <- as.character(Hmat[,1])
        Hmat <- Hmat[,-1]
    } else {
        signatureNames <- as.character(1:nrow(Hmat))
    }
    rownames(Hmat) <- signatureNames
    
    # Calculate frequencies
    #Hmat <- t(t(Hmat) / (colSums(Hmat)))
    
    currentSignature <-  signifCor[i,"signature"]
    
    signatureText <- paste0(basename(Hfile), "_signature", currentSignature)
    
    
    # Reading age and foucusing on individuals for which we have age info
    age <- readMetadataGTEX("AGE")[sraToGtex(colnames(Hmat)), , drop = F ]
    colnames(Hmat) <- sraToGtex(colnames(Hmat))
    age <- age[rownames(age) %in% colnames(Hmat),,drop=F]
    Hmat <- Hmat[, rownames(age)]
    nInds <- ncol(Hmat)
    
    age_all <- c(age_all, age[,1])
    value <- c(value, as.numeric(Hmat[currentSignature, ]))
    signature <- c(signature, rep(signatureText, nInds))
    signifCor[i, "signature"] <- signatureText 
    
}

forScatter <- data.frame(signature = signature, value = value, age = age_all)
##signifCor$label <- paste0("pval = ", signif(signifCor$pvalue,3), "\n",
#                          "cor = ", signif(signifCor$cor,3), "\n"
#                          )

#signifCor$x <- signifCor$y <- Inf

p <- scatter(forScatter, x = "age", y = "value", facet_x = "signature", nrowFactor = ceiling(nrow(signifCor)/3), regression = T) + theme_noGrid()
#p <- ggplot(forScatter, aes(x = age, y = value)) +
#    geom_point() +
#    stat_smooth(method = "lm") + 
#    geom_text(aes(x = x, y = y, label = label), data = signifCor, hjust = 1, vjust = 1) + 
#    facet_wrap(~signature, scales = "free", ncol = 3) +
#    theme_bw() 

#p_nrow <- ceiling(nrow(signifCor) / 3)
#p_ncol <- ifelse(nrow(signifCor) < 3, nrow(signifCor), 3)

p_nrow =  ceiling(nrow(signifCor)/3)
p_ncol = 3
    
ggsave(outScatterFile, p, height = p_nrow * PANEL_HEIGHT, width = p_ncol * PANEL_WIDTH)
ggsave(outHistFile, pHist)


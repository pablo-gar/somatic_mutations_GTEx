# Takes all the mutation counts for a given tissue/MAF combination
# the read count table and the metatadata table.
#
# Does a linear regression gets the pvalues and coeffiecients of all the metavariables
# This is to know how much the mutation counts are influenced by external factors

# Usage
# Rscript mutationVsReads.R outScatter.pdf readCountTable.txt count1.txt [count2.txt] [...] 
   
if(!dir.exists("../../R"))
   stop("Can't find the share R path, make sure to run script from directory where it lives")

source("../../R/plots.R")
source("../../R/ggthemes.R")
source("../../R/gtex.R", chdir = T)

args <- commandArgs(T)
outScatter <- args[1]
readCountFile <- args[2]
transDiversityFile <- args[3]
args <- args[4:length(args)]

#readCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"
#transDiversityFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/transcriptome_shannon_diversity.txt"
#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/Whole_Blood/n6_0.0_0.7"
#args <- file.path(workingDir, list.files(workingDir))

metadataCols <- c("SMRIN", #RIN
                  "SMTSISCH", #Ischemic time
                  "SMTSPAX" #Time in fixative
                  )

indInfo <- c("AGE", "GENDER", "BMI", "RACE")

factorVariables <- c("GENDER")
featureNames <- c("RIN", "Ischemic time", "Time in fixative", "Age", "Self-defined race", "Gender", "BMI", "Uniquely mapped reads", "Transcriptome diversity", "Genotype PC1", "Genotype PC2", "Genotype PC3")
names(featureNames) <- c("SMRIN", "SMTSISCH", "SMTSPAX", "AGE", "RACE", "GENDER", "BMI", "n_uniqueMapped", "transcriptome_diversity", "C1", "C2", "C3")

artifacts <- c("SMRIN", "SMTSISCH", "SMTSPAX",  "n_uniqueMapped", "transcriptome_diversity")
noArtifacts <- names(featureNames)[!names(featureNames) %in% artifacts]

#------------------------------------#
# METHODS

# Read count table should be produced with ../bin/getReadCounts
readGeneCounts <- function(x) {
    readCountTable <- read.table(x, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    return(readCountTable)
}

readTransDiversity <- function(x) { 
    transDiversity <- read.table(x, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    return(transDiversity)
}

# Read mutation count tables
readMutations <- function(files) {
    countTable <- data.frame(sample = "", mutations = rep(0, length(args)), stringsAsFactors = F)

    for (i in 1:length(files)) {
        mutation <- read.table(files[i], sep = "\t", stringsAsFactors = F, header = F)
        countTable[i, "sample"] <- gsub(".txt", "", basename(files[i]))
        countTable[i, "mutations"] <- sum(mutation)
    }
    
    return(countTable)
    
}


#------------------------------------#
# MAIN

# Read mutation count tables
countTable <- readMutations(args)
sraIds <- countTable$sample
gtexIds <- sraToGtex(countTable$sample, formatOut = "short")
gtexIdsLong <- sraToGtex(countTable$sample, formatOut = "long")

rownames(countTable) <-  gtexIdsLong

# Read gene counts file
readCountTable <- readGeneCounts(readCountFile)

# Read transcriptome diversity
transDiversity <- readTransDiversity(transDiversityFile)

# Read metadata
sampleMeta <- readSampleAnnotationGTEX(metadataCols)
indMeta <- readMetadataGTEX(indInfo)
genotypePCA <- t(readGenotypePCA())


# Getting ids for which we have info for all
allInfo <- sraIds %in% rownames(readCountTable) & 
    gtexIdsLong %in% rownames(sampleMeta) & 
    gtexIds %in% rownames(indMeta) & 
    gtexIds %in% rownames(genotypePCA) &
    gtexIdsLong %in% rownames(transDiversity)

sraIds <- sraIds[allInfo]
gtexIds<- gtexIds[allInfo]
gtexIdsLong <- gtexIdsLong[allInfo]

# Merging tables
countTable <- countTable[gtexIdsLong,]

countTable <- cbind(countTable, readCountTable[sraIds,"n_uniqueMapped", drop = F])
countTable <- cbind(countTable, transDiversity[gtexIdsLong, , drop = F])
countTable <- cbind(countTable, sampleMeta[gtexIdsLong,])
countTable <- cbind(countTable, indMeta[gtexIds,])
countTable <- cbind(countTable, genotypePCA[gtexIds,1:3])
countTable <- countTable[,-1]

# Eliminate columns with all na
countTable <- countTable[,colSums(is.na(countTable)) != nrow(countTable)]

# Doing lm and getting pvalues for each coefficient

if(nrow(countTable) > 1) {
    pvalues <- -log10(summary(lm(mutations ~ ., data = countTable))[["coefficients"]][,"Pr(>|t|)"])
    pvalues <- data.frame(log10_p = pvalues, feature = featureNames[names(pvalues)])
    pvalues <- pvalues[-1,] 
    
    p <- ggplot(pvalues, aes(x = feature, y = log10_p)) + geom_bar(stat = "identity")
    p <- p + xlab("feature") + ylab("-log10(pvalue)") + coord_flip() + theme_grid_x()
    
    
    # Writing table
    rownames(pvalues) <- pvalues$feature
    pvalues <- pvalues[,"log10_p", drop = F]
    pvalues <- t(pvalues)
    write.table(pvalues, paste0(outScatter, ".pvals.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

} else {
    p <- ggplot(data.frame(x = 1, y = 1), aes(x = x, y = y)) + geom_point()
}

ggsave(outScatter, p, width = 5, height = 5)


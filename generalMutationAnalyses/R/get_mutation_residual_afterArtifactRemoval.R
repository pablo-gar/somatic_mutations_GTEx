# Takes all the mutation counts for a given tissue/MAF combination
# the read count table and the metatadata table.
#
# Gets:
#   1. Table with mutations counts and covariates
#   2. Table with mutations resitudals after removal of artifacts, it appends the other covariates too
# Gets a table with the residuals from a linear regression using artifacts as coefficients
# 
# it also appends non-artifact columns 

# Usage
# Rscript get_mutation_residual_afterArtifactRemoval.R outTable.txt readCountTable.txt count1.txt [count2.txt] [...] 
   
if(!dir.exists("../../R"))
   stop("Can't find the share R path, make sure to run script from directory where it lives")

source("../../R/plots.R")
source("../../R/ggthemes.R")
source("../../R/gtex.R", chdir = T)

library("tidyr")
library("dplyr")

args <- commandArgs(T)
outTable <- args[1]
outTable_original <- args[2]
readCountFile <- args[3]
transDiversityFile <- args[4]
args <- args[5:length(args)]

#readCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"
#transDiversityFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/transcriptome_shannon_diversity.txt"
#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/Whole_Blood/n6_0.0_0.5"
#args <- file.path(workingDir, list.files(workingDir))


allMuts <- c("C_A", "C_G", "C_T", "T_A", "T_C", "T_G", "all")
metadataCols <- c("SMRIN", #RIN
                  "SMTSISCH", #Ischemic time
                  "SMTSPAX" #Time in fixative
                  )

indInfo <- c("AGE", "RACE", "GENDER", "BMI")

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
    countTable <- list()

    for (i in 1:length(files)) {
        mutation <- read.table(files[i], sep = "\t", stringsAsFactors = F, header = F)
        allCount <- sum(mutation)
        
        # Narrow format
        colnames(mutation) <- c("A", "T", "G", "C")
        mutation$from <- c("T", "C")
        mutation <- gather(mutation, "to", "mutations", -from)
        
        # Expand to single row each mutation one column
        mutation$mut <- paste0(mutation$from, "_", mutation$to)
        mutation <- 
            mutation %>%
            select(mut, mutations) %>%
            spread(mut, mutations)
        
        mutation$all <- allCount
        
        # Append sample name
        mutation$sample <- gsub(".txt", "", basename(files[i]))
        
        countTable[[i]] <- mutation
    }
    
    countTable <- do.call(rbind, countTable)
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

# Eliminate columns with all na
countTable <- countTable[,colSums(is.na(countTable)) != nrow(countTable)]

# Writes original counts
write.table(countTable, outTable_original, sep = "\t", quote = F, col.names = T, row.names = F)

if(nrow(countTable) > 1) {
    
    for(mut in allMuts) {
        form <- as.formula(paste(mut, "~ ."))
        artResiduals <- resid(lm(form, data = countTable[,colnames(countTable) %in% c(mut, artifacts)]))
        countTable[names(artResiduals), mut] <- artResiduals
    }    
        
} else {
    p <- ggplot(data.frame(x = 1, y = 1), aes(x = x, y = y)) + geom_point()
}

# Write table
countTable <- countTable[ ,!colnames(countTable) %in% artifacts]
write.table(countTable, outTable, sep = "\t", quote = F, col.names = T, row.names = F)

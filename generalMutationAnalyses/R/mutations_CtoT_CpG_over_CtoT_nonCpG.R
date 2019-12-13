# Calculates the ratio of C>T mutations in CpG vs C>T in non CpGs
# Takes a context matrix and ONLY WORKS WITH CONTEXT LENGTHS OF 3

# Usage
# Rscript mutationVsExpression.R mapFile1.txt [mapFile2.txt] [...]

library("dplyr")
args <- commandArgs(T)
contextFile <- args[1]
outFile <- args[2]


#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.5/"
#mutationFiles <- file.path(workingDir, list.files(workingDir))

#----------------------------------#
# MAIN
#----------------------------------#
allMutations <-  read.table(contextFile, sep = "\t", stringsAsFactors = F, colClasses = c("character", "character", "numeric"))

# Select C>T
C_to_T <- allMutations[allMutations[,1] == "11",]
C_to_T <- C_to_T[!grepl("N", C_to_T[,2]),]
CpG_rows <- grepl("..G", C_to_T[,2])
results <- data.frame(CpG = sum(C_to_T[CpG_rows, 3]), non_CpG = sum(C_to_T[!CpG_rows, 3]))
results$fc_CpG_nonCpG <- log2(results$CpG / results$non_CpG)

# Writing ouputs
write.table(results, outFile, sep = "\t", quote = F, row.names = F, col.names = T)

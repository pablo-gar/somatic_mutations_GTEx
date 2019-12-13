# Creates a histogram of distances between mutations and also a table of quantiles of the same values
# Usage
# Rscript mutationVsExpression.R mapFile1.txt [mapFile2.txt] [...]

library("dplyr")
args <- commandArgs(T)
mutationFile <- args[1]
outFile <- args[2]


#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.5/"
#mutationFiles <- file.path(workingDir, list.files(workingDir))

#----------------------------------#
# MAIN
#----------------------------------#
allMutations <-  read.table(mutationFile, sep = "\t", stringsAsFactors = F)
allMutations$distance <- 0

# Calculating distances
for(chr in unique(allMutations[,1])) {
    
    currentChr <- allMutations[ allMutations[,1] == chr, ] 
    currentChr <- currentChr[order(currentChr[,2]), ]
    currentChr[-nrow(currentChr), "distance"] <- currentChr[-1, 2] - currentChr[-nrow(currentChr), 2]
    
    allMutations[ allMutations[,1] == chr, ] <- currentChr
}

# CALCULATE AVERAGE AND VAR AND SD
distanceStats <-
    allMutations %>%
    summarise(median = median(distance), mean = mean(distance), sd = sd(distance), var = var(distance), total = n())

# Writing ouputs
write.table(distanceStats, outFile, sep = "\t", quote = F, row.names = F, col.names = T)

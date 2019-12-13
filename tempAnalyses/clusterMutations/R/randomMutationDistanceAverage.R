# Takes a mutation file that has at least 2 columns, chrom and position.
# Select n random rows and calculates median, mean, sd, and var of distance between them

# Usage
# Rscript mutationVsExpression.R mapFile1.txt [mapFile2.txt] [...]

library("dplyr")
args <- commandArgs(T)
mutationFile <- args[1]
nRandomMuts <- as.integer(args[2])
permutations <- as.integer(args[3])
outFile <- args[4]


#----------------------------------#
# MAIN
#----------------------------------#
allMutations <-  read.table(mutationFile, sep = "\t", stringsAsFactors = F)
allDistanceStats <- list()
for (i in 1:permutations) {
    randomMutations <- allMutations[sample(nrow(allMutations), nRandomMuts), ]
    randomMutations$distance <- 0

    # Calculating distances
    for(chr in unique(randomMutations[,1])) {
        
        currentChr <- randomMutations[ randomMutations[,1] == chr, ] 
        currentChr <- currentChr[order(currentChr[,2]), ]
        currentChr[-nrow(currentChr), "distance"] <- currentChr[-1, 2] - currentChr[-nrow(currentChr), 2]
        
        randomMutations[ randomMutations[,1] == chr, ] <- currentChr
    }

    # CALCULATE AVERAGE AND VAR AND SD
    distanceStats <-
        randomMutations %>%
        summarise(median = median(distance), mean = mean(distance), sd = sd(distance), var = var(distance))
    
    allDistanceStats[[i]] <- distanceStats
}

# Writing ouputs
allDistanceStats <- do.call(rbind, allDistanceStats)
allDistanceStats$n <- nRandomMuts
write.table(allDistanceStats, outFile, sep = "\t", quote = F, row.names = F, col.names = T)

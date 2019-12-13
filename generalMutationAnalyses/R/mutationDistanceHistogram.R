# Creates a histogram of distances between mutations and also a table of quantiles of the same values
# Usage
# Rscript mutationVsExpression.R mapFile1.txt [mapFile2.txt] [...]

library("ggplot2")
source("../../R/ggthemes.R", chdir = T)
args <- commandArgs(T)
outPrefix <- args[1]
mutationFiles <- args[2:length(args)]


#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.5/"
#mutationFiles <- file.path(workingDir, list.files(workingDir))

#----------------------------------#
# METHODS
#----------------------------------#
readMutations <- function(mutationFiles) {
    mutationsAll <- list()
    for(i in 1:length(mutationFiles)) {
        
        mutationFile <- mutationFiles[i]
        
        currentSample <- gsub(".txt", "", basename(mutationFile))
        mutations <- read.table(mutationFile, sep = "\t", stringsAsFactors = F)
        mutations$sample <- currentSample
        
        mutationsAll[[i]] <- mutations
        
    }

    mutationsAll <- do.call(rbind, mutationsAll)
    
    return(mutationsAll)
}

#----------------------------------#
# MAIN
#----------------------------------#
allMutations <- readMutations(mutationFiles)
allMutations$distance <- 0

# Calculating distances
for(currentSample in unique(allMutations$sample)) {
    for(chr in unique(allMutations[,1])) {
        
        currentChr <- allMutations[ allMutations[,1] == chr & allMutations$sample == currentSample,] 
        currentChr <- currentChr[order(currentChr[,2]), ]
        currentChr[-nrow(currentChr), "distance"] <- currentChr[-1, 2] - currentChr[-nrow(currentChr), 2]
        
        allMutations[ allMutations[,1] == chr & allMutations$sample == currentSample,] <- currentChr
    }
}


# Getting plot of only up to 85 percentile
allMutations <- allMutations[allMutations$distance > 0,]
maxLimit <- quantile(allMutations$distance, 0.85)
toPlot <-  allMutations[allMutations$distance < maxLimit, ]


p <- ggplot(toPlot, aes(x = distance)) +
geom_histogram() +
theme_grid_y() +
theme(axis.text.x = element_text(angle = 30, hjust = 1))


# Writing ouputs

ggsave(paste0(outPrefix, "distHist.pdf"))

stats <- summary(allMutations$distance)
write.table(data.frame(stat = names(stats), value = as.vector(stats)), paste0(outPrefix, "stats.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

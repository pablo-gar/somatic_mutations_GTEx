# Makes map of median distances between adjacent basepairs after having
# create many mutation maps by excluding different type of clustered mutations
# based on the size of the cluster and the number of mutations per cluster

library("purrr")
library("tidyr")
library("dplyr")
library("ggplot2")
source("../../../R/plots.R") 

main <- function(cmdArgs = commandArgs(T)) {
    
    tissueFolder <- "/scratch/users/paedugar/somaticMutationsProject/clusterCombinations/distanceStats/Whole_Blood/n6_0.0_0.5"
    tissueFolderExpected <- "/scratch/users/paedugar/somaticMutationsProject/clusterCombinations/distanceStatsExpected/Whole_Blood/n6_0.0_0.5"
    sampleFiles <-  list.files(tissueFolder)
    
    distance <- readDistanceFiles(sampleFiles, tissueFolder, tissueFolderExpected)
    distance <- distance[!is.na(distance$expectedMedian),]
    distance <- mutate(distance, fold_change = log2(median/expectedMedian))
    

    fc_stats <-
        distance %>%
        group_by(length, nMut) %>%
        summarise(median_mutations = median(total), mean_fc = mean(fold_change), sd_fc = sd(fold_change)) %>%
        mutate(ymax = mean_fc + sd_fc, ymin = mean_fc - sd_fc) %>%
        ungroup()
    
    ggplot(fc_stats, aes(x = nMut)) +
        geom_pointrange(aes(y = mean_fc, ymin = ymin, ymax = ymax, colour = median_mutations)) +
        geom_text(aes(y = mean_fc + 2, label = median_mutations), colour = "grey60", vjust = 1) +
        facet_grid(length~.) + 
        theme_bw() +
        theme(legend.position = "top")

    
}

readDistanceFiles <- function(sampleFiles, dirDist, dirDistExpected) {
    
    # File notation is [sample]_length=[number]_nMuts=[number].txt
    
    distance <- map_dfr(sampleFiles, function(x) {
                            out <-  read.table(file.path(dirDist, x), header = T, stringsAsFactors = F, sep = "\t")
                            out$sample <- gsub("(\\w+?)_.+", "\\1", x)
                            out$length <- as.integer(gsub(".+length.(\\d+).+", "\\1", x))
                            out$nMut <- as.integer(gsub(".+nMuts.(\\d+).+", "\\1", x))
                            
                            # Getting expected
                            if(file.exists(file.path(dirDistExpected, x))) {
                                expected <-  read.table(file.path(dirDistExpected, x), header = T, stringsAsFactors = F, sep = "\t")
                                expected <- mean(expected$median)
                            } else {
                                expected <- NA
                            }
                            
                            out$expectedMedian <- expected
                            
                            return(out)
                           }
                        )
    
    return(distance)

}

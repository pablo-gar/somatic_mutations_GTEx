# Takes files containing distances in each line and plots a histogram of them
# Usage
# Rscript plot_exonDistance.R outplot.pdf /scratch/users/paedugar/somaticMutationsProject/mutationCount/closestExon/Pancreas/n6_0.0_0.5/*

library("tidyr")
library("dplyr")
library("purrr")
library("ggplot2")

cmdArgs <- commandArgs(T)
outPlot <- cmdArgs[1]
distanceFiles <- cmdArgs[2:length(cmdArgs)]

# Plot preferences
pointColour <- "grey50" 
xmin <- 0
xmax <- 100
psize <- 0.2

width <- 8.5
height <- 2.5

#folderDistances <- "/scratch/users/paedugar/somaticMutationsProject/exonBoundaryDistance/closestExon_merged/Adipose_Subcutaneous/n6_0.0_0.7"
#folderDistances <- "/scratch/users/paedugar/somaticMutationsProject/exonBoundaryDistance/closestExon_merged/Whole_Blood_EXO/n6_0.0_0.7"
#folderDistances <- "/scratch/users/paedugar/somaticMutationsProject/exonBoundaryDistance/closestExon_merged/Lung/n6_0.0_0.7"
#outPlot <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/distanceClosestExonHist/Whole_Blood_EXO-n6_0.0_0.7_distanceHist.pdf"
#distanceFiles <- file.path(folderDistances, list.files(folderDistances))

#---------------------------
# METHODS

readDistances <- function(x) {
    
    result <- map_dfr(x, function(y) { 
                       distances <- read.table(y, sep = "\t", header = F, stringsAsFactors = F)
                       distances[2,2:3] <- distances[2,2:3] + distances[1,2:3]
                       distances <- distances[distances[,1] >= 0, ]
                       colnames(distances) <- c("distance", "count_all", "count_mutations")
                       distances$mut_rate <- (distances$count_mutations + 0.0001) / distances$count_all
                       distances$sample <- basename(y)
                       distances <- distances[distances$count_all > 0, ]
                       return(distances)
                      }
                    )
    return(result)
}

# Takes a vector of values, returns the difference of the k elemnts uptream minus the k elements dowstream (including n)
differenceGroup <- function(x, k) {
    
    result <- rep(0, length(x))
    for(i in 1:length(result)) {
        
        strUps <- ifelse( i + 1 > length(x), length(x), i + 1)
        endUps <- ifelse( i + k > length(x), length(x), i + k)
        
        srtDws <- ifelse(i - k + 1 < 1, 1, i - k + 1)
        endDws <- i
        
        result[i] <- mean(x[srtDws:endDws]) - mean(x[strUps:endUps])
    }
    
    return(result)
}



#---------------------------
# MAIN
distances <- readDistances(distanceFiles)

# Eliminate crazy outliers
meanMutRate <- tapply(distances$mut_rate, distances$sample, mean)
outliersQuant <-  quantile(meanMutRate,  c(0.01, 0.99))
outliers <- names(meanMutRate[meanMutRate <= outliersQuant[1] | meanMutRate >= outliersQuant[2]])
distances <- distances[!distances$sample %in% outliers,]


# Get mean counts, window values, and difference accross windows
smoothWindow <- 4
distanceCountMean <- 
    distances %>%
    group_by(distance) %>%
    summarise(mean_n = mean(mut_rate), sd_n = sd(mut_rate), 
              total = n()
              ) %>%
    ungroup() %>%
    mutate(#Window values
           diff_n = (differenceGroup(mean_n, k = smoothWindow)),
           # Difference to upstream window
           diff_n_next = c(diff_n[2:length(diff_n)] - diff_n[1:(length(diff_n) - 1)], 0),
           )
                                    



cutoff <- distanceCountMean$distance[which.min(distanceCountMean$diff_n_next[1:25])]
toPlot <- distanceCountMean[distanceCountMean$distance >= xmin & distanceCountMean$distance <= xmax,]

p <- ggplot(toPlot, aes(x = distance, y = mean_n)) +
        annotate("rect", xmin = cutoff - smoothWindow + 0.5, xmax = cutoff + 0.5, ymin = min(distanceCountMean$mean_n), ymax = Inf, fill = "#f2aa82", alpha = 0.4) + 
        annotate("rect", xmax = cutoff + smoothWindow + 0.5, xmin = cutoff + 0.5, ymin = min(distanceCountMean$mean_n), ymax = Inf, fill = "#448cff", alpha = 0.4) + 
        annotate("text", x = cutoff + 1.5, y = min(distanceCountMean$mean_n), label = paste("d =", cutoff), vjust = 0, hjust = 0) +
        geom_vline(xintercept = cutoff + 0.5, linetype = "dotted") +
        geom_pointrange(aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n), colour =  pointColour, size = psize) +
        geom_line(colour = "grey20") + 
        xlim(xmin, xmax) + 
        ylab("Mean mutations") +
        xlab("Bps to closest exon boundary") +
        theme_bw()

# FOR EXOME
#p <- p + ylim(0, 0.0035)

ggsave(outPlot, p, width = width, height = height)


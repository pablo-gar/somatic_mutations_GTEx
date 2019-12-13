# Compares C to T ratios in CpG context vs not by excluding different types of clusters


library("purrr")
library("tidyr")
library("dplyr")
library("ggplot2")
source("../../../R/plots.R") 

main <- function(cmdArgs = commandArgs(T)) {
    
    tissueFolder <- "/scratch/users/paedugar/somaticMutationsProject/clusterCombinations/C_to_T_ratios/Whole_Blood/n6_0.0_0.5"
    sampleFiles <-  list.files(tissueFolder)
    
    c_to_t <- readFiles(sampleFiles, tissueFolder)
    c_to_t <- mutate(c_to_t, fold_change = log2(CpG/non_CpG))
    

    fc_stats <-
        c_to_t %>%
        group_by(length, nMut) %>%
        summarise(mean_fc = mean(fold_change), sd_fc = sd(fold_change)) %>%
        mutate(ymax = mean_fc + sd_fc, ymin = mean_fc - sd_fc) %>%
        ungroup()
    
    ggplot(fc_stats, aes(x = nMut)) +
        geom_pointrange(aes(y = mean_fc, ymin = ymin, ymax = ymax)) +
        facet_grid(length~.) + 
        ylab("Mean log2 fold-change C>T (CpG / non-CpG") +
        xlab("Number of mutations in clusters") +
        theme_bw() +
        theme(legend.position = "top")

    
}

readFiles <- function(sampleFiles, dirDist) {
    
    # File notation is [sample]_length=[number]_nMuts=[number].txt
    
    c_to_t <- map_dfr(sampleFiles, function(x) {
                            out <-  read.table(file.path(dirDist, x), header = T, stringsAsFactors = F, sep = "\t")
                            out$sample <- gsub("(\\w+?)_.+", "\\1", x)
                            out$length <- as.integer(gsub(".+length.(\\d+).+", "\\1", x))
                            out$nMut <- as.integer(gsub(".+nMuts.(\\d+).+", "\\1", x))
                            
                            return(out)
                           }
                        )
    
    return(c_to_t)

}

#' Takes one or more map files with a filter column, with filters
#' separated by ";" and computes basic stats
#'
#' @param 0-based column number for the filter column
#' @param filter names seprated by ",". Order matters
#' @param one or more map files separated by spaces
#' @example
#' Rscript filter_stats.R 7 blacklisted_region,rna_edit,splicing_junction_error,sequencing_error,clustered_mutation out.pdf /scratch/users/paedugar/somaticMutationsProject/mutationCount/filtered/Adipose_Subcutaneous/n6_0.0_0.5/*


source("../../R/misc.R")
source("../../R/ggthemes.R")
library("purrr")
library("dplyr")
library("tidyr")
library("ggplot2")

main <- function(cmdArgs = commandArgs(T)) {
    
    column_filter <- as.integer(cmdArgs[1]) + 1 
    filter_names <- parseArg(cmdArgs[2], sep = ",")
    outTable <- cmdArgs[3]
    outPlot <- cmdArgs[4]
    map_files <- cmdArgs[5:length(cmdArgs)]
    
    
    #column_filter <- as.integer(7) + 1 
    #filter_names <- parseArg("blacklisted_region,rna_edit,splicing_junction_error,sequencing_error,bcf_read_position_bias,bcf_mapping_quality_bias,bcf_base_quality_bias,bcf_mapping_quality_vs_strand_bias,bcf_variant_distance_bias,clustered_mutation", sep = ",")
    #outTable <- "out.txt"
    #outPlot <- "out.pdf"
    #map_files <- list.files("/scratch/users/paedugar/somaticMutationsProject/mutationCount/filtered/Lung/n6_0.0_0.7", full = T)
   
   
    # Read filter counts per file
    filter_numbers <- read_map_files(map_files, column_filter, filter_names)
    filter_numbers <- filter_numbers[, colnames(filter_numbers) != "clustered_mutation"]
    
    # Sum across samples
    grand_total <-  filter_numbers[,"total"]
    for(i in 3:ncol(filter_numbers)) {
        filter_numbers[,i] <- (filter_numbers[,i-1]- filter_numbers[,i])
    }
    
    filter_numbers[,2:ncol(filter_numbers)]  <- filter_numbers[,2:ncol(filter_numbers)] / grand_total * 100
    # Calculate means per filter
    filter_numbers_means <- c(mean(grand_total), colMeans(filter_numbers[,2:ncol(filter_numbers)]))
    filter_numbers_means <- matrix(filter_numbers_means, nrow = 1)
    colnames(filter_numbers_means) <- c("Total", colnames(filter_numbers)[2:ncol(filter_numbers)])
    
    
    
    
    filter_stats <-
        filter_numbers %>%
        gather(key = "filter", value = "mutations", -sample) %>%
        group_by(filter)  %>%
        summarise(mean_left = mean(mutations), mean_eliminated = 100 - mean_left, lower_conf = bootstrap_confidence_interval(mutations, base::mean, 10e3)[1], upper_conf = bootstrap_confidence_interval(mutations, base::mean, 10e3)[2]) %>%
        ungroup() %>%
        gather("mean_type", "mean", mean_left, mean_eliminated) 
    
    filter_stats$filter <- factor(filter_stats$filter, levels = rev(c("total", filter_names)), ordered = T)
    # Only takes unique info for extra added details
    filter_stats_info <-  filter_stats[filter_stats$mean_type == "mean_left",]
    
    p <- ggplot(filter_stats, aes( x = filter, y = mean)) +
        geom_bar(aes(fill = mean_type), stat = "identity", width = 0.8) +
        geom_errorbar(aes(ymin = lower_conf, ymax = upper_conf), data = filter_stats_info, width = 0.2) + 
        geom_text(aes(label = round(mean,2), y = mean - 2), colour = "#4299ff", data = filter_stats_info, hjust = 1) +
        geom_text(aes(label = filter), y = 0, colour = "#4299ff", fontface = "bold", data = filter_stats_info, hjust = 0) +
        ylab("Mutations left (%)") + xlab("Cummulative filter substraction") + 
        scale_fill_manual(values = c("#e5f1ff", "#b5d9ff")) + 
        theme_void() + 
        theme(legend.position = "none") +
        coord_flip()
    
    ggsave(outPlot, p)
    write.table(filter_numbers_means, outTable, sep = "\t", col.names = T, row.names =F, quote = F)
    
}

read_map_files <- function(map_files, column_filter, filter_names) {
    
    result <- map_dfr(map_files, column_filter = column_filter, filter_names = filter_names, 
                      .f = function(x, column_filter, filter_names){
                          
                          current_map <- read.table(x, sep = "\t", stringsAsFactors = F)
                          current_map <- current_map[,column_filter]
                          
                          filter_counts <- list()
                          filter_counts[["sample"]] <- basename(x)
                          filter_counts[["total"]] <- length(current_map)
                          for(c_filter in filter_names) {
                              filter_l <- grepl(c_filter, current_map)
                              filter_counts[[c_filter]] <- sum(filter_l)
                              current_map <- current_map[!filter_l]
                          }
                          filter_counts = do.call(data.frame, filter_counts)
                          return (filter_counts)
                      })
    
    return(result)
                                                               
}

main()

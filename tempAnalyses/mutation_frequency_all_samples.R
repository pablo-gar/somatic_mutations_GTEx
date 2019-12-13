# Performs a histogram of mutation frequency across entire population
# Input is the output of
# /home/users/paedugar/scripts/FraserLab/somaticMutationsProject/generalMutationAnalyses/bin/print_mut_freq2 /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/*/n6_0.0_0.7/* > /scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/mutation_freq_table.txt
#
# /home/users/paedugar/scripts/FraserLab/somaticMutationsProject/generalMutationAnalyses/bin/print_mut_freq2 /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/*/n6_0.0_0.7/* > /scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/mutation_freq_table.txt
# Rscript mutation_frequency_all_samples.R /scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/ /scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/mutation_freq_table.txt
#
# NO EXO
# mv /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Whole_Blood_EXO/ /scratch/users/paedugar/somaticMutationsProject/mutationCount/Whole_Blood_EXO_map
# /home/users/paedugar/scripts/FraserLab/somaticMutationsProject/generalMutationAnalyses/bin/print_mut_freq2 /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/*/n6_0.0_0.7/* > /scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/noExo/mutation_freq_table.txt
# mv /scratch/users/paedugar/somaticMutationsProject/mutationCount/Whole_Blood_EXO_map/ /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Whole_Blood_EXO 
#
# Rscript mutation_frequency_all_samples.R /scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/noExo/ /scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/noExo/mutation_freq_table.txt


library("ggplot2")
source("../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    
    out_prefix <- cmdArgs[1]
    freq_table_file <- cmdArgs[2]
    
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/"
    #freq_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/mutation_freq_table.txt"
    
    # Read
    freq_table <- read.table(freq_table_file, sep = "\t", stringsAsFactors = F, header = F)
    freq_table$V1 <- 100 * freq_table$V1
    freqs <- rep(freq_table$V1, freq_table$V3)
    
    # Do histogram
    p <- plot_histogram(freqs, bins = 100, xlabel = "Mutation frequency across samples (%) -- all mutations")
    p2 <- plot_histogram(freq_table$V1, bins = 100, xlabel = "Mutation frequency across samples (%) -- unique mutations")
    
    # Plot total vs unique freqs
    if (max(freq_table$V1) > 5) {
        p3 <- plot_total_unique(freqs, freq_table$V1, values = unique(c(seq(0,1, 0.1), seq(1,5, 0.2), seq(5, max(freq_table$V1), 1))))
        ggsave(paste0(out_prefix, "unique_vs_total.pdf"), p3, width = 18)
    }
    
    ggsave(paste0(out_prefix, "histogram.pdf"), p)
    ggsave(paste0(out_prefix, "unique_mutations_histogram.pdf"), p2)
    
    # Save stats
    stats <- data.frame(stat = c("first_quartile", "second_quartile", "mean_freq", "median_freq", "unique_mutations", "total_mutations", "total_samples", "max", "min"),
                        value = c(quantile(freqs, 0.25), (quantile(freqs, 0.75)), mean(freqs), median(freqs), nrow(freq_table), length(freqs), unique(freq_table$V2), max(freq_table$V1), min(freq_table$V1))
                        )
    
    total_muts <- data.frame(uniq = nrow(freq_table), total = sum(freq_table[,3]))
        
    write.table(stats, paste0(out_prefix, "stats.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
    write.table(total_muts, paste0(out_prefix, "stats_total_muts.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
    
    
}

plot_total_unique <- function(a, b, values, xlab = "Total mutations freq", ylab = "Unique mutations freq") {
    toPlot <- data.frame(a = 0, b = 0, value = values)
    length_a <- length(a)
    length_b <- length(b)
    toPlot$a <- sapply(values, function(x) sum(a < x) / length_a)
    toPlot$b <- sapply(values, function(x) sum(b < x) / length_b)
    
    ggplot(toPlot, aes(x = a, y = b)) +
        geom_point() +
        geom_line() +
        geom_text(aes(label = value), vjust = 0, size = 3) +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw()
    
    
}


plot_histogram <- function(x, bins = 30, xlabel) {
    
    numbers <- paste0("average = ", signif(mean(x), 2), "\n",
                      "median = ", signif(median(x), 2), "\n",
                      "1st quartile = ", signif(quantile(x, 0.25), 2), "\n",
                      "3d quartile = ", signif(quantile(x, 0.75), 2), "\n",
                      "n = ", length(x), "\n"
                      )
                      
    
    p <- ggplot(data.frame(x = x), aes(x = x)) +
            geom_histogram(bins = bins, colour = "black", fill = NA) +
            annotate("text", x = Inf, y = Inf, label = numbers, hjust = 1, vjust = 1) +
            xlab(xlabel) +
            theme_noGrid()
        
    return(p)
}


main()

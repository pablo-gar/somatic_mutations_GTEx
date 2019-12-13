# takes output of print_mut_freq and plots percentages of mutation plots divided by
# pop frequency bins. It also produces a histogram of variant allele frequency (pop level)
# of all mutations
#
# out of print_mut_freq:
# 1. frequency of mutation in population
# 2. chr
# 3. pos
# 4. ref
# 5. alt

library("dplyr")
library("ggplot2")

main <- function(cmdArgs = commandArgs(T)) {
    
    #ARGS
    vaf_cutoff <- as.numeric(cmdArgs[1])
    mut_table <- cmdArgs[2]
    outPlot <- cmdArgs[3]
    outPlot_histogram <- cmdArgs[4]
    
    #mut_table <-"/scratch/users/paedugar/somaticMutationsProject/mutationCount/pop_freq/Brain_Cerebellum-n6_0.0_0.7.txt"
    
    #VARS
    ranges <- c(0, 0.05, 0.1, 0.2, 0.5, 1)
    rangesNames <- paste(ranges[-length(ranges)], "-", ranges[-1])
    
    freq_table <- read.table(mut_table, sep = "\t", stringsAsFactors = F, header = F)
    colnames(freq_table) <- c("pop_frequency", "chrom", "pos", "from", "to")
    
    #Plots histograms
    p_hist <- ggplot(freq_table, aes(x = pop_frequency)) +
        geom_histogram(fill = NA, colour = "grey30") +
        xlab("VAF in population") +
        geom_vline(xintercept = vaf_cutoff, linetype = "dashed", colour = "grey50") +
        theme_bw()
    
    
    # Read and group mutations
    freq_table$mut <- getMutations(paste0(freq_table$from, ">", freq_table$to), sep = ">")
    freq_table$group  <- findInterval(freq_table$pop_frequency, ranges, rightmost.closed = T)
    freq_table$group  <- findInterval(freq_table$pop_frequency, ranges, rightmost.closed = T)
    freq_table$group <- rangesNames[freq_table$group]
    
    # Calculate frequencies per bin
    freq_types <- 
        freq_table %>%
        group_by(mut, group) %>%
        summarise(counts = n()) %>%
        ungroup() %>%
        group_by(group) %>%
        mutate(freq = counts / sum(counts)) %>%
        ungroup() %>%
        mutate(group = factor(group, levels = rangesNames, order = T))
    
    p <- ggplot(freq_types, aes(x = group, y = freq, fill = mut)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("#3BC9F3", "#2B2E34", "#FC3218", "#CAD0CE", "#9CD169", "#F1CAC9")) +
            theme_bw()
        
    ggsave(outPlot, p)
    ggsave(outPlot_histogram, p_hist)
    
    
}

getMutations <- function(x, sep = ">") {
    
    bases <- setNames(c("A", "T", "G", "C"), c("T", "A", "C", "G"))
    muts <- data.frame(from = rep(c("A", "T", "G", "C"), each = 4), to = c("A", "T", "G", "C"), stringsAsFactors = F)
    to_trans <- muts$from %in% c("A", "G")
    muts_names <- muts
    muts_names$to[to_trans] <-bases[muts_names$to[to_trans]]
    muts_names$from[to_trans] <-bases[muts_names$from[to_trans]]
    
    mutKey <- setNames(paste0(muts_names$from, sep, muts_names$to), paste0(muts$from, sep, muts$to))
    
    return(mutKey[x])
}

main()

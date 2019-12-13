library("ggplot2")
library("dplyr")
library("tidyr")
main <- function(cmdArgs = commandArgs(T)) {
    
    #mutation_dir <- "/scratch/users/paedugar/somaticMutationsProject/chromatin/mutation_exon_signal/Spleen/n6_0.0_0.7/"
    #mutation_files <- file.path(mutation_dir, list.files(mutation_dir)[1:50])
    #chromatin_file <- "/scratch/users/paedugar/somaticMutationsProject/chromatin/roadmap_exon_signal_compiled/E106.txt"
    
    outPlot <- cmdArgs[1]
    outTable <- cmdArgs[2]
    chromatin_file <- cmdArgs[3]
    mutation_files <- cmdArgs[4:length(cmdArgs)]
    
    # Reads chromatin files
    chromatin <- read.table(chromatin_file , sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
    mods <- colnames(chromatin)
    chromatin$exon_length <- 0
    
    # Reads mutations files
    for(i in mutation_files) {
        mutations <- read.table(i, sep = "\t", stringsAsFactors = F, header = T)
        
        # Add mutations
        chromatin$mut <- 0
        chromatin[mutations$exon, "mut"] <- mutations$mutations
        colnames(chromatin)[ncol(chromatin)] <- paste0("mut_", basename(i))
        
        # Add bps sequenced
        chromatin$bps_seq <- 0
        chromatin[mutations$exon, "bps_seq"] <- mutations$total_bps_seq
        colnames(chromatin)[ncol(chromatin)] <- paste0("bps_seq_", basename(i))
        
        chromatin[mutations$exon, "exon_length"] <- mutations$exon_length
    }
    
    # Get exons for which we have info in all marks
    all_marks <- rowSums(chromatin[,mods] != 0) == length(mods)
    
    # Calculate sum of mutations and bps sequenced
    chromatin_compiled <- chromatin[all_marks,c("exon_length", mods)]
    chromatin_compiled$mutations <- rowSums(chromatin[all_marks, grep("mut", colnames(chromatin))])
    chromatin_compiled$total_bps_seq <- rowSums(chromatin[all_marks, grep("bps_seq", colnames(chromatin))])
    
    chromatin_compiled <- chromatin_compiled[chromatin_compiled$mutations != 0,]
    chromatin_compiled$mutations <- log2(chromatin_compiled$mutations)
    chromatin_compiled$total_bps_seq <- log2(chromatin_compiled$total_bps_seq)
    windows <-  round(nrow(chromatin_compiled) / (nrow(chromatin_compiled) * 0.005))
    chromatin_compiled <- get_residuals(chromatin_compiled, windows, T)
    
    # Get residuals median per window
    
    # Do linear regression and correlations
    #chromatin_compiled$mut_rate <- log(chromatin_compiled$mutations / chromatin_compiled$total_bps_seq)
    chromatin_compiled$mut_rate <- chromatin_compiled$mutations - chromatin_compiled$median_mut_group
    chromatin_compiled <- chromatin_compiled[,c("mut_rate", mods)]
    coefficents <- summary(lm(mut_rate ~ . , data = chromatin_compiled))$coefficients
    cors <- cor(chromatin_compiled)[mods, "mut_rate"]
    coefficents <- as.data.frame(cbind(coefficents, cor = c(0,cors)))
    coefficents$mod <- rownames(coefficents)
    
    # Plot
    toPlot <-
        chromatin_compiled %>%
        gather("mod",  "fc", -mut_rate) %>%
        group_by(mod) %>%
        mutate(quartile = factor(findInterval(fc, quantile(fc, seq(0, 1, length.out = 5)), all.inside = T))) %>%
        ungroup()
    
    p <- ggplot(toPlot, aes(x = quartile, y = mut_rate)) +
            #geom_boxplot() +
            geom_violin()+
            #facet_grid( ~ mod, scales = "free") +
            facet_wrap( ~ mod, ncol = 3, scales = "free") +
            theme_bw()
        
    ggsave(outPlot, p)
    write.table(coefficents, outTable, sep = "\t", row.names = F, quote = F)
    
    
}

get_residuals <- function(chromatin_compiled, windows = 100, by_quantile = T) {
    
    if(by_quantile) {
        groups <- quantile(chromatin_compiled$total_bps_seq, seq(0,1, length.out=windows))
    } else {
        groups <- seq(min(chromatin_compiled$total_bps_seq), max(chromatin_compiled$total_bps_seq), length.out = windows)
    }
    chromatin_compiled$group <- findInterval(chromatin_compiled$total_bps_seq, groups, all.inside = T)
    
    medians <-
        chromatin_compiled %>%
        group_by(group) %>%
        summarise(median_mut = median(mutations)) %>%
        ungroup() %>%
        as.data.frame()
    rownames(medians) <- medians$group
    
    chromatin_compiled$median_mut_group <- medians[as.character(chromatin_compiled$group), "median_mut"]
    
    return(chromatin_compiled)
    
}
main()

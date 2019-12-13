# Plots percentages of mutations that are also observed in the cosmic database
library("ggplot2")
library("dplyr")
source("../../../R/ggthemes.R")
source("../../../R/misc.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    mutationDir <- cmdArgs[1]
    maf <- cmdArgs[2]
    outplot_file <- cmdArgs[3]
    
    HEIGHT  <- 4
    WIDTH <- 12
    
    #mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/map"
    #annotated_mut_file <- "/scratch/users/paedugar/somaticMutationsProject/cancer/dndsout_tissuesMerged_all/1_all_mutations_annotated_dnds.txt.gz"
    #maf <- "n6_0.0_0.7"
    #outplot_file <-  "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/plots/n6_0.0_0.7/vaf_comparison.pdf"
    
    
    
    # Read annotated mutations
    annotated_mutations <- read.table(annotated_mut_file, header = T, sep = "\t", stringsAsFactors = F)
    annotated_mutations <- annotated_mutations[,c("chr", "pos", "ref", "mut", "impact", "gene_name")]
    annotated_mutations <- unique(annotated_mutations)
    annotated_mutations$id <- paste0("chr", annotated_mutations$chr, ".", annotated_mutations$pos, ".", annotated_mutations$ref, ".", annotated_mutations$mut)
    annotated_mutations <- annotated_mutations[annotated_mutations$impact %in% c("Synonymous", "Missense", "Nonsense"),]
    
    # Read mutations and select those that are annotated
    mutations <- read_cosmic_maps(mutationDir, maf)
    mutations$id <- paste0(mutations$chr, ".", mutations$start, ".", mutations$ref, ".", mutations$alt) 
    
    mutations_an <- inner_join(mutations, annotated_mutations, by = "id")
    mutations_an$cancerous <- ifelse(mutations_an$alt_cosmic == mutations_an$alt, "Seen_in_cancer", "Not_in_cancer")
    mutations_an$VAF <- mutations_an$alt_count / mutations_an$coverage
    
    # Getting median values
    mutation_stats <- mutations_an %>%
        group_by(impact, cancerous)  %>% 
        summarise(Median_VAF = median(VAF), lo_conf = bootstrap_confidence_interval(VAF, median, 1000)[1], up_conf = bootstrap_confidence_interval(VAF, median, 1000)[2]) %>%
        ungroup()
    
    mutation_stats$impact <- factor(mutation_stats$impact, levels =  c("Synonymous", "Missense", "Nonsense"), ordered = T)
    mutation_stats$cancerous <- factor(mutation_stats$cancerous, levels =  c("Seen_in_cancer", "Not_in_cancer"), ordered = T)
    
    # Getting stats 
    mutations_signif <- mutations_an %>%
        group_by(impact) %>%
        summarise(p = signif(wilcox.test(VAF ~ cancerous)[["p.value"]], 2), p_label = ifelse(p == 0, "p < 1e-300", paste0("p = ", p))) %>%
        ungroup()
        
    p <- ggplot(mutation_stats, aes(x = impact, y = Median_VAF)) +
        geom_bar(aes(fill = cancerous), colour = "black", stat = "identity", position = position_dodge(width = 0.9)) +
        geom_errorbar(aes(ymax = up_conf, ymin = lo_conf, fill = cancerous), position = position_dodge(width = 0.9), width = 0.3) + 
        xlab("") + ylab("Median VAF") +
        geom_text(aes(label = p_label), y = Inf, data = mutations_signif, vjust = 1) + 
        theme_grid_y() + 
        theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1))
    
    ggsave(outplot_file, p, width = 3, height = 5)
    
}

read_cosmic_maps <- function (mutationDir, maf) {
    
    results <- list()
    i <- 1 
    for(tissue in list.dirs(mutationDir, recursive = F)) {
        if (grepl("EXO", tissue))
            next
        for(current_sample in list.files(file.path(tissue, maf))) {
            flush.console()
            print(paste(tissue, current_sample))
            info <- file.info(file.path(tissue, maf, current_sample))
            if(any(info$size == 0))
                next
            
            current <- read.table(file.path(tissue, maf, current_sample), sep = "\t", stringsAsFactors = F, header = F)
            colnames(current) <- c("chr_cosmic", "start_cosmic", "end_cosmic", "ref_cosmic", "alt_cosmic", "strand", "chr", "start", 
                                   "end", "ref", "alt", "context", "coverage", "alt_count")
            current$tissue <- basename(tissue)
            current$sample <- gsub(".txt", "", basename(current_sample))
            results[[i]] <- current
            i <- i + 1
        }
    }
    
    results <- do.call(rbind, results)
    
    return(results)
    
}

main()

# Takes cancer mutations and returns a table with mutations that are likely false positives
library("ggplot2")
library("dplyr")
library("tidyr")
source("../../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    cutoff_hist <- as.numeric(cmdArgs[1])
    mutation_dir <- cmdArgs[2]
    outplot_hist <- cmdArgs[3]
    out_table <- cmdArgs[4]
    
    
    # DEBUG
    #mutation_dir <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverMutations/mutations/n6_0.0_0.7"
    #cutoff_hist <- 3.5
    #out_table <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverMutations/mutations_false_positives/n6_0.0_0.7/cancer_mutations_false_positives.txt"
    #outplot_hist <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverMutations/mutations_false_positives/n6_0.0_0.7/cancer_mutations_false_positives_histogram.pdf"
    
    # Reading files
    mutations <- read_mutation_files(mutation_dir)
    mutations <- mutations[!grepl("EXO", mutations$tissue),]
    
    # Plot special histogram, for each mutation calculate a value of [# times mutation is observed] / [# unique mutations in associated gene]
    # it returns those indices for mutations that fall above a defined cutoff
    mutation_freq <- plot_freq_mutations(mutations, cutoff = cutoff_hist)
    
    mutations_to_eliminate <- mutations[mutation_freq$indices_to_eliminate,]
    mutations_to_eliminate <- mutations_to_eliminate[,c("chr", "pos", "ref", "mut", "gene_name")]
    mutations_to_eliminate <- unique(mutations_to_eliminate) 
    
    p_eliminated_mutations_hist <- mutation_freq$p
    
    ggsave(outplot_hist, p_eliminated_mutations_hist)
    write.table(mutations_to_eliminate, out_table, sep = "\t", quote = F, col.names = F, row.names = F)
    
    
    # Second method
    mutation_freq <- plot_scatter_mutations(mutations, cutoff = 15)
    
    mutations_to_eliminate <- mutations[mutation_freq$indices_to_eliminate,]
    mutations_to_eliminate <- mutations_to_eliminate[,c("chr", "pos", "ref", "mut", "gene_name")]
    mutations_to_eliminate <- unique(mutations_to_eliminate) 
    
    p_eliminated_mutations_hist <- mutation_freq$p
    
    ggsave(paste0(outplot_hist, ".old.pdf"), p_eliminated_mutations_hist)
    write.table(mutations_to_eliminate, paste0(out_table, ".old"), sep = "\t", quote = F, col.names = F, row.names = F)
    
    
}

read_mutation_files <- function(x) {
    
    x <- file.path(x, list.files(x))
    results <-list()
    for(i in seq_along(x)) {
        current <- read.table(x[i], sep = "\t", stringsAsFactors = F, header = T)
        current$tissue <- gsub(".txt", "", basename(x[i]))
        results[[i]] <- current
    }
    results <- do.call(rbind, results)
    return(results)
}

# Plot special histogram, for each mutation calculate a value of [# times mutation is observed] / [# unique mutations in associated gene]
# it returns those indices for mutations that fall above a defined cutoff
plot_freq_mutations <- function(mutations, cutoff = 0) {
    
    mutations$mut <- paste0(mutations$chr, ".", mutations$pos, ".", mutations$ref, ".", mutations$mut, ".", mutations$gene_name)
    
    # count unique mutations per gene
    muts_per_gene <- mutations%>%
        group_by(gene_name) %>%
        summarise(total = length(unique(mut))) %>%
        ungroup() %>%
        as.data.frame()
    rownames(muts_per_gene) <- muts_per_gene$gene_name
    
    # for each unique mutation count occurrances
    mut_count <- mutations %>%
        group_by(mut) %>%
        summarise(total = n()) %>%
        ungroup() %>%
        as.data.frame()
    rownames(mut_count) <- mut_count$mut
    
    # Across all mutations append their value of count of mutation / unique counts in corresponding gene
    mutations$counts_over_unique <- 0
    for(i in 1:nrow(mutations)) {
        total_mut_count <-  mut_count[mutations[i, "mut"], "total"]
        unique_muts_in_gene <- muts_per_gene[mutations[i, "gene_name"], "total"]
        mutations[i, "counts_over_unique"] <- total_mut_count / unique_muts_in_gene
    }
    
    
    p <- ggplot(mutations, aes(x=counts_over_unique)) +
            geom_histogram(fill = NA, colour = "black", stat = "count") +
            xlab("Mutations in cancer driver genes") +
            geom_vline(xintercept = cutoff, linetype = "dashed", colour = "red") +
            theme_grid_y()
        
    indices_to_eliminate = mutations$counts_over_unique >= cutoff
    return(list(p = p, indices_to_eliminate = indices_to_eliminate))
    
}

plot_scatter_mutations <- function(mutations, cutoff = 0) {
    
    mutations$mut <- paste0(mutations$chr, ".", mutations$pos, ".", mutations$ref, ".", mutations$mut, ".", mutations$gene_name)
    
    toPlot <- tapply(mutations$mut, mutations$gene_name, function(x) {
                          aaa <- as.data.frame(table(x))
                          aaa$n = log10(length(x))
                          aaa$Freq = aaa$Freq/length(x)
                          aaa  
                         })
    
    toPlot <- do.call(rbind, toPlot)
    
    toPlot$to_eliminate <- F
    toPlot[toPlot$n > 1 & toPlot$Freq > .30, "to_eliminate"] <- T
    toPlot[toPlot$n > 3 & toPlot$Freq > .10, "to_eliminate"] <- T
    
    indices_to_eliminate <- mutations$mut %in% toPlot[toPlot$to_eliminate, 1]
    
    
    p <- ggplot(toPlot, aes(x=n, y = Freq)) +
        geom_point(aes(colour = to_eliminate)) +
        annotate("text", label = paste0("keeping = ", sum(!indices_to_eliminate), "\nremoving = ", sum(indices_to_eliminate)), x = Inf, y = Inf, hjust = 1, vjust = 1) + 
        theme_bw()
    
        
    return(list(p = p, indices_to_eliminate= indices_to_eliminate))
    
}

main()

# Plots the results of mutationsDriverMutationsOverlap.R after
# being applied to all tissues.
# Makes counts per tissue, top genes mutated in each tissue and 
# ratio of syn vs nonSyn mutations
library("ggplot2")
library("gplots")
library("dplyr")
library("tidyr")
library("ggrepel")
source("../../../R/ggthemes.R")
source("../../../R/misc.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    out_prefix <- cmdArgs[1]
    dnds_file <- cmdArgs[2]
    mutation_dir <- cmdArgs[3]
    expression_dir <- cmdArgs[4]
    expected_mutations_file <- cmdArgs[5]
    
    
    # DEBUG
    #dnds_file <-  "/scratch/users/paedugar/somaticMutationsProject/cancer/dndsout_tissuesMerged_noFalsePos_long_list/n6_0.0_0.7/dndsout_sel_loc.txt"
    #mutation_dir <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverMutations_long_list/mutations/n6_0.0_0.7"
    #expression_dir <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverMutations_long_list/expression/n6_0.0_0.7"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverGeneralPlots_long_list/n6_0.0_0.7/"
    #expected_mutations_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/observed_predicted_counts_per_tissue/predicted_observed.txt"
    #oncogenic_muts_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/oncoKB/allAnnotatedVariants_pointMut_oncogenic.txt"
    #gene_driver_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/tcga_driverCancerGenes.txt"
    
    # Reading files
    #oncogenic_muts <- read.table(oncogenic_muts_file, sep = "\t", header = T, stringsAsFactors = F)
    #rownames(oncogenic_muts) <- paste0(oncogenic_muts$Gene, ".", oncogenic_muts$Alteration)
    
    mutations <- read_mutation_files(mutation_dir)
    mutations <- mutations[!grepl("EXO", mutations$tissue),]
    
    # Appending oncogene vs TSG status
    driver_status <- read_driver_status(gene_driver_file)
    mutations$driver_status <-  driver_status[mutations$gene_name]
    
    
    expected_mutations <- read.table(expected_mutations_file, sep = "\t", stringsAsFactors = F, header = T)
    expected_mutations <- expected_mutations[expected_mutations$mut_type == "mutations",]
    rownames(expected_mutations) <- expected_mutations$tissue
    
    expression_genes <- read_mutation_files(expression_dir)
    expression_genes <- expression_genes[!grepl("EXO", expression_genes$tissue),]
    
    dnds <- read.table(dnds_file, sep = "\t", header = T, stringsAsFactors = F)
    
    
    # Plot special histogram, for each mutation calculate a value of [# times mutation is observed] / [# unique mutations in associated gene]
    # it returns those indices for mutations that fall above a defined cutoff
    mutation_scatter <- plot_scatter_mutations(mutations, cutoff = 15)
    #mutations <- mutations[mutation_scatter$indices,]
    
    mutation_freq <- plot_freq_mutations(mutations, cutoff = 3.5)
    mutations <- mutations[mutation_freq$indices,]
    
    p_eliminated_mutations_hist <- mutation_freq$p
    p_eliminated_mutations_scatter <- mutation_scatter$p
    
    
    # PLOTS
    oncokb_table <- make_oncokb(mutations) 

    
    #Bar plots of mutations on cancer driver genes per tissue
    p_genes_per_tissue <- plot_nGene_per_tissue(mutations, expected_mutations, F)
    p_genes_per_tissue_mut_rate <- plot_nGene_per_tissue(mutations, expected_mutations, T)
    
    p_scatter_genes <- plot_top_genes_per_tissue(mutations)
    
    p_impacts <- plot_genes_mutation_impacts(mutations)
    p_expression_all_tissues <- plot_gene_expression_all_tissues(expression_genes, mutations)
    p_dnds <- plot_dnds(dnds, mutations, n_mutations_cutoff = 5)
    
    p_dnds_box <- plot_dnds_boxplots(dnds, mutations, n_mutations_cutoff = 5, out_prefix = out_prefix, n_perm = 10000, FUN_stat = mean) 
    p_dnds_box_oncogene <- plot_dnds_boxplots(dnds, mutations[mutations$driver_status== "oncogene",], n_mutations_cutoff = 5, out_prefix = out_prefix, n_perm = 10000, FUN_stat = mean) 
    p_dnds_box_tsg <- plot_dnds_boxplots(dnds, mutations[mutations$driver_status== "tsg",], n_mutations_cutoff = 5, out_prefix = out_prefix, n_perm = 10000, FUN_stat = mean) 
    
    
    # Get number of mutations
    pdf(paste0(out_prefix, "presence_heatmap.pdf"), width = 11, height = 11)
    plot_gene_presence_heatmap(mutations, expected_mutations, F)
    dev.off()
    
    # Do tissue-wise mutation rate per gene
    pdf(paste0(out_prefix, "presence_heatmap_mutation_rate_perTissue"), width = 11, height = 11)
    plot_gene_presence_heatmap(mutations, expected_mutations, T)
    dev.off()
    
    ggsave(paste0(out_prefix, "histogram_mutation_counts_over_unique_mut_in_gene.pdf"), p_eliminated_mutations_hist, width = 13, height = 4)
    ggsave(paste0(out_prefix, "scatter_mutation_counts_over_unique_mut_in_gene.pdf"), p_eliminated_mutations_scatter, width = 13, height = 4)
    ggsave(paste0(out_prefix, "driverCount_per_tissue.pdf"), p_genes_per_tissue, width = 13, height = 4)
    ggsave(paste0(out_prefix, "driverCount_per_tissue_mut_rate.pdf"), p_genes_per_tissue_mut_rate, width = 13, height = 4)
    ggsave(paste0(out_prefix, "scatterDriver_per_tissue.pdf"), p_scatter_genes, width = 18, height = 18)
    ggsave(paste0(out_prefix, "impact_mutations.pdf"), p_impacts, width = 4, height = 8)
    ggsave(paste0(out_prefix, "gene_expression_per_gene.pdf"), p_expression_all_tissues, width = 2.5, height = 8)
    ggsave(paste0(out_prefix, "dnds_mutations.pdf"), p_dnds, width = 3, height = 10)
    ggsave(paste0(out_prefix, "dnds_mutations_boxplot.pdf"), p_dnds_box, width = 3, height = 7)
    ggsave(paste0(out_prefix, "dnds_mutations_boxplot_oncogne.pdf"), p_dnds_box_oncogene, width = 3, height = 7)
    ggsave(paste0(out_prefix, "dnds_mutations_boxplot_tsg.pdf"), p_dnds_box_tsg, width = 3, height = 7)
    
    # Write mutations and also write the tables for the cBioportal lolliplots 
    out_dir <- file.path(out_prefix, "mutation_tables")
    dir.create(out_dir)
    write.table(mutations, paste0(out_prefix, "all_working_mutations.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(oncokb_table, paste0(out_prefix, "oncokb_mutations.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
    for(gene in unique(mutations$gene_name)) {
        current_table <- mutations[mutations$gene_name == gene, c("sample", "tissue", "chr", "pos", "pos", "ref", "mut")]
        colnames(current_table) <- c("Sample_ID", "Cancer_Type", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Variant_Allele")
        write.table(current_table, file.path(out_dir, paste0(gene, "_.txt")), sep = "\t", quote = F, row.names = F, col.names = T)
    }
           
    
}

read_driver_status <- function(x) {
    
    x <- read.delim(x, sep = "\t", stringsAsFactors = F, header = T)
    x$Tumor_suppressor_or_oncogene_prediction <- gsub(".+_", "", x$Tumor_suppressor_or_oncogene_prediction)
    
    status <- x %>%
        group_by(Gene, Tumor_suppressor_or_oncogene_prediction) %>%
        summarise(count = n()) %>%
        ungroup()
    
    status <- status[order(status$Gene, -status$count),]
    status <- status[!duplicated(status$Gene),]
    status <- setNames(status[, 2, drop = T], status[, 1, drop = T])
    return(status)
        
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
            scale_x_discrete(breaks=seq(0,max(mutations$counts_over_unique), by = 2)) +
            theme_grid_y()
        
    indices_to_keep = mutations$counts_over_unique < cutoff
    return(list(p = p, indices = indices_to_keep))
    
}

# Plot special scatter, freq mutation vs total number of mutations in that gene
plot_scatter_mutations <- function(mutations, cutoff = 0) {
    
    mutations$mut <- paste0(mutations$chr, ".", mutations$pos, ".", mutations$ref, ".", mutations$mut, ".", mutations$gene_name)
    
    toPlot <- tapply(mutations$mut, mutations$gene_name, function(x) {
                          aaa <- as.data.frame(table(x))
                          aaa$n = log10(length(x))
                          aaa$Freq = aaa$Freq/length(x)
                          aaa  
                         })
    
    toPlot <- do.call(rbind, toPlot)
    toPlot$gene <- gsub("(\\w+)\\.", "\\1", rownames(toPlot))
    
    toPlot$to_eliminate <- F
    toPlot[toPlot$n > 1 & toPlot$Freq > .30, "to_eliminate"] <- T
    toPlot[toPlot$n > 3 & toPlot$Freq > .10, "to_eliminate"] <- T
    
    indices_to_keep <- !mutations$mut %in% toPlot[toPlot$to_eliminate, 1]
    
    
    p <- ggplot(toPlot, aes(x=n, y = Freq)) +
        geom_point(aes(colour = to_eliminate)) +
        annotate("text", label = paste0("keeping = ", sum(indices_to_keep), "\nremoving = ", sum(!indices_to_keep)), x = Inf, y = Inf, hjust = 1, vjust = 1) + 
        theme_bw()
    
        
    return(list(p = p, indices = indices_to_keep))
    
}
    
plot_nGene_per_tissue <- function(mutations, expected_mutations, normalize_expected = F, expected_mutations_column = "observed") {

    gene_per_tissue <- mutations %>%
        group_by(tissue) %>%
        #summarise(n_genes = length(unique(gene_name))) %>%
        summarise(n_genes = length(gene_name)) %>%
        ungroup() 
    
    if(normalize_expected) {
        gene_per_tissue$n_genes <- gene_per_tissue$n_genes / expected_mutations[gene_per_tissue$tissue, expected_mutations_column]
        ylab = '[Mutations in driver genes] / [Total mutations in tissue expected by seq depth]'
    } else {
        ylab = 'Mutations in driver genes'
    }
    
    # Get significant tissues
    signif_tissues <- tissue_driverMut_significant(mutations, expected_mutations, expected_mutations_column)
    signif_tissues$label <- labelPvalues(signif_tissues$fdr)
    
    gene_per_tissue <- mutate(gene_per_tissue, tissue = factor(tissue, levels = tissue[order(-n_genes)], ordered = T))
    gene_per_tissue$fdr_label <- signif_tissues[as.character(gene_per_tissue$tissue), "label"]
    
    p <- ggplot(gene_per_tissue, aes(x = tissue, y = n_genes)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = fdr_label), hjust = 0.5, vjust = 0) +
            geom_hline(aes(yintercept = mean(n_genes)), linetype = "dashed") +
            ylab(ylab) +
            theme_grid_y() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
            #coord_flip()
        
    return(p) 
}

#' Gets tissues that are significanltly more or less mutated (in their driver genes) than the averages
tissue_driverMut_significant <- function(mutations, expected_mutations, expected_mutations_column) {
    
    total_observed <- sum(expected_mutations[,expected_mutations_column])
    total_cancer_muts <- nrow(mutations)
    rate <- total_cancer_muts / total_observed
    
    mutations_per_tissue <- tapply(mutations[,1], mutations$tissue, length)
    
    # Calculating pvalue
    expected_mutations$driver_muts <- mutations_per_tissue[rownames(expected_mutations)]
    expected_mutations$pvalue_lowtail <- pbinom(q = expected_mutations$driver_muts, size = expected_mutations[,expected_mutations_column], prob = rate, log.p = T)
    expected_mutations$pvalue_hightail <- pbinom(q = expected_mutations$driver_muts, size = expected_mutations[,expected_mutations_column], prob = rate, log.p = T, lower.tail = F)
    
    # Calculating fdr
    expected_mutations$fdr_lowtail <- p.adjust(10^expected_mutations$pvalue_lowtail, method = "BH")
    expected_mutations$fdr_hightail <- p.adjust(10^expected_mutations$pvalue_hightail, method = "BH")
    
    # Get lowest fdr
    expected_mutations$fdr <- ifelse(expected_mutations$fdr_lowtail < expected_mutations$fdr_hightail, expected_mutations$fdr_lowtail, expected_mutations$fdr_hightail)
    
    return(expected_mutations)
    
}

plot_top_genes_per_tissue <- function(mutations) {
    
    driver_mutations <- mutations %>%
        group_by(tissue, gene_name) %>%
        summarise(n_muts = n()) %>%
        ungroup() %>%
        group_by(tissue) %>%
        mutate(gene_ord = rank(n_muts, ties.method="first"), gene_ord_rev = rank(-n_muts, ties.method="first")) %>%
        ungroup()
    
    driver_mutations$gene_label <- driver_mutations$gene_name
    driver_mutations$gene_label[driver_mutations$gene_ord_rev >6] <- ""
    
    gene_per_tissue <- mutations %>%
        group_by(tissue) %>%
        summarise(n_genes = length(unique(gene_name))) %>%
        ungroup() 
    tissues <- gene_per_tissue$tissue[order(-gene_per_tissue$n_genes)]
    
    driver_mutations$tissue <- factor(driver_mutations$tissue, levels = tissues, ordered = T)
    
    
    p <- ggplot(driver_mutations, aes(x = gene_ord, y = n_muts)) +
            facet_wrap(~tissue, scales = "free") +
            geom_point() +
            geom_text(aes(label = gene_label), hjust = 1) +
            #geom_text_repel(aes(label = gene_label), direction = "y") +
            ylab("Mutations") +
            theme_classic() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    return(p)
    
}

plot_genes_mutation_impacts <- function(mutations, dnds) {
    
    impact_cols <- c("#df7373", "#8D5FD3", "#f9d0a7", "#9DAC93")
    
    mutation_impacts <- mutations %>%
        group_by(gene_name, impact) %>%
        summarise(counts = n()) %>%
        ungroup() %>%
        group_by(gene_name) %>%
        mutate(percent = counts/sum(counts) * 100) %>%
        ungroup()
    
    # Getting gene counts
    gene_count <- mutations %>%
        group_by(gene_name) %>%
        summarise(mut_count = n()) %>%
        ungroup()
    
    # Get gene order
    genes_ordered <- get_gene_orders(mutations)
    
    mutation_impacts$gene_name <- factor(mutation_impacts$gene_name, levels = rev(genes_ordered), ordered = T)
    
    # Plot percentage of different mutation impacts
    p <- ggplot(mutation_impacts, aes(x = gene_name, y = percent)) +
            geom_bar(aes(fill = impact), stat = "identity") +
            scale_fill_manual(values = impact_cols) +
            geom_text(aes(label = paste("n =", mut_count)), fontface = "bold.italic", colour = "white", y = 95, hjust = 1, data = gene_count) +
            ylab("Mutations (%)") +
            xlab("") +
            theme_classic() +
            coord_flip() +
            theme(legend.position = "top")
        
    return(p)
}

plot_gene_expression_all_tissues <- function(expression_genes, mutations) {
    
    genes_ordered <- get_gene_orders(mutations)
    
    toPlot <- expression_genes %>%
        filter(gene_name %in% mutations$gene_name) %>%
        group_by(gene_name) %>%
        summarise(value = log2(mean(total_reads)), x = "a") %>%
        ungroup()
    
    toPlot$gene_name <- factor(toPlot$gene_name, levels = rev(genes_ordered), ordered = T)
    
    p <- ggplot(toPlot, aes(x = x, y = gene_name)) + 
            geom_tile(aes(fill = value)) +
            scale_fill_continuous(name = "No. reads") +
            theme_bw()
        
    return(p)
    
}

get_dnds_mean <- function(dnds, get_log2 = T, n_mutations_cutoff = 5) {
    
    # Getting wmin average and non sense
    dnds_signif_mis <- get_dnds_good_genes(dnds, "wmis_loc", n_mutations_cutoff = n_mutations_cutoff)$wmis_loc
    dnds_signif_non <- get_dnds_good_genes(dnds, "wnon_loc", n_mutations_cutoff = n_mutations_cutoff)$wnon_loc
    
    if(get_log2) {
        dnds_signif_mis <- log2(dnds_signif_mis)
        dnds_signif_non <- log2(dnds_signif_non)
    }
    
    dnds_mis_conf <- bootstrap_confidence_interval(dnds_signif_mis, mean)
    dnds_non_conf <- bootstrap_confidence_interval(dnds_signif_non, mean)
    
    dnds_mis_conf_median <- bootstrap_confidence_interval(dnds_signif_mis, median)
    dnds_non_conf_median <- bootstrap_confidence_interval(dnds_signif_non, median)
    
    dnds_mis_mean <- data.frame(mean = mean(dnds_signif_mis), conf_lo = dnds_mis_conf[1], conf_up = dnds_mis_conf[2], 
                                median = median(dnds_signif_mis), conf_lo_median = dnds_mis_conf_median[1], conf_up_median = dnds_mis_conf_median[2]
                                )
    dnds_non_mean <- data.frame(mean = mean(dnds_signif_non), conf_lo = dnds_non_conf[1], conf_up = dnds_non_conf[2], 
                                median = median(dnds_signif_non), conf_lo_median = dnds_non_conf_median[1], conf_up_median = dnds_non_conf_median[2]
                                )
    
    return(list( mis = dnds_mis_mean, non = dnds_non_mean))
    
}


# Test how many time a randomly selected group of genes is above the median dnds
get_pvalue_dnds_median <- function(dnds, genes, measure = "wmis_loc", n_perm = 10e5, lower.tail = F, n_mutations_cutoff = 5, get_log2 = T, FUN_stat = median) {
    
    
    # Real median
    if(get_log2)
        dnds[,measure] <- log2(dnds[,measure])
    
    real_m <- do.call(FUN_stat, list(x = dnds[dnds$gene_name %in% genes, measure]))
    
    dnds <- get_dnds_good_genes(dnds, measure, n_mutations_cutoff = n_mutations_cutoff)
    # Permutated medians
    permuted_m <- sapply(rep(length(genes), n_perm), dnds = dnds, measure = measure,
                         FUN = function(x, dnds, measure) do.call(FUN_stat, list(x=(dnds[sample(nrow(dnds), x, replace = F), measure])))
                         )
    
    if(lower.tail) {
        pvalue <- sum(permuted_m <= real_m) / n_perm
    } else {
        pvalue <- sum(permuted_m >= real_m) / n_perm
    }
    
    return(pvalue)
    
}

# Get genes from the dnds table for which we feel confident
get_dnds_good_genes <- function(dnds, measure = "wmis_loc", n_mutations_cutoff = 5) {
    
    if(measure == "wmis_loc") {
        mutColumn <- "n_mis"
        #dnds[,measure] <- dnds[p.adjust(dnds$pmis_loc, method = "BH") < 0.01, measure] 
    } else {
        mutColumn <- "n_non"
    }
    
    # only mutations with at least one syn and one from measure of interest
    dnds <- dnds[dnds$n_syn > 0 & dnds[,mutColumn] > 0, ] 
    
    # Only genes that have the number of mutations cutoff
    dnds <- dnds[dnds$n_syn +  dnds[,mutColumn] >= n_mutations_cutoff, ]
    
    # if measure equals 0 lets create the lowest available
    dnds[dnds[,measure] == 0, measure] <- 1 / (max(dnds[,measure]) + 1)
    
    return(dnds)
    
}

plot_dnds_boxplots <- function(dnds, mutations, n_mutations_cutoff = 5, out_prefix = out_prefix, FUN_stat = median, n_perm = 10000) {
    # get dnds mean and conf intervals
    dnds_mean <- get_dnds_mean(dnds, n_mutations_cutoff = n_mutations_cutoff)
    dnds_mis_mean <-  dnds_mean[[1]]
    dnds_non_mean <-  dnds_mean[[2]]
    
    # get only genes of interest 
    dnds_all <- dnds
    dnds <- dnds[dnds$gene_name %in%  mutations$gene_name,]
    
    # only genes with more than n mutations
    gene_count <- table(mutations$gene_name)
    gene_count_cutoff <- gene_count[gene_count >= n_mutations_cutoff]
    dnds <- dnds[dnds$gene_name %in%  rownames(gene_count_cutoff),]
    
    toPlot <- 
        dnds %>%
        gather(type, value, wmis_loc, wnon_loc) %>%
        mutate(value = log2(value))
    
    # Gather wmis if we have more than 1 missense and wnon if we have more than 1 nonsense, and for all at least 1 syn mutation
    toPlot <- toPlot[! (toPlot$type == "wmis_loc" & toPlot$n_mis == 0),]
    toPlot <- toPlot[! (toPlot$type == "wnon_loc" & toPlot$n_non == 0),]
    toPlot <- toPlot[toPlot$n_syn >0, ]
    
    # Getting pvalues
    wmis_genes <- unique(toPlot[toPlot$type == "wmis_loc", "gene_name"])
    wnon_genes <- unique(toPlot[toPlot$type == "wnon_loc", "gene_name"])
    pvalues <- data.frame(type = c("wmis_loc", "wnon_loc"),
                          label = c(get_pvalue_dnds_median(dnds_all, wmis_genes, measure = "wmis_loc", n_mutations_cutoff = n_mutations_cutoff, FUN_stat = FUN_stat, n_perm = n_perm),
                                   get_pvalue_dnds_median(dnds_all, wnon_genes, measure = "wnon_loc", n_mutations_cutoff = n_mutations_cutoff, FUN_stat = FUN_stat, n_perm = n_perm)
                                   ),
                          value = Inf
                          )
    pvalues$label <- paste("p =", signif(pvalues$label, 1))
    
    write.table(toPlot, paste0(out_prefix, "dnds_boxplot_values.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(rbind(dnds_mis_mean, dnds_non_mean), paste0(out_prefix, "dnds_boxplot_values_genome_wide.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
    
    ggplot(toPlot, aes(x = type, y = value)) +
        geom_boxplot(colour = c("#8AC4FF", "#FFC39E"), fill = NA, lwd = 1.5) +
        annotate("rect", colour = "#8AC4FF", fill = "#8AC4FF", ymin = dnds_mis_mean$conf_lo, ymax = dnds_mis_mean$conf_up, xmin = -Inf, xmax = Inf) +
        annotate("rect", colour = "#FFC39E", fill = "#FFC39E", ymin = dnds_non_mean$conf_lo, ymax = dnds_non_mean$conf_up, xmin = -Inf, xmax = Inf) +
        geom_hline(yintercept = dnds_mis_mean$mean, colour = "#0073e6", linetype = "dashed") +
        geom_hline(yintercept = dnds_non_mean$mean, colour = "#e65800", linetype = "dashed") +
        stat_summary(fun.y=mean, colour = c("#8AC4FF", "#FFC39E"), geom="point", shape=18, size=5, show_guide = FALSE) +
        geom_text_repel(aes(label = gene_name), size = 3.2) +
        geom_text(aes(label = label), data = pvalues, vjust = 1) +
        ylab("log2(dN/dS)") + 
        theme_grid_y()
    
    
}

plot_dnds <- function(dnds, mutations, n_mutations_cutoff = 4) {
    
    genes_ordered <- get_gene_orders(mutations)
    
    
    # get dnds mean and conf intervals
    dnds_mean <- get_dnds_mean(dnds, n_mutations_cutoff = n_mutations_cutoff)
    dnds_mis_mean <-  dnds_mean[[1]]
    dnds_non_mean <-  dnds_mean[[2]]
    
    # Ordering data for plot
    dnds <- dnds[dnds$gene_name %in%  genes_ordered,]
    dnds$gene_name <- factor(dnds$gene_name, levels = rev(genes_ordered), ordered = T)
    
    toPlot <- 
        dnds %>%
        gather(type, value, wmis_loc, wnon_loc) %>%
        mutate(value = log2(value))
        #mutate(value = (value))
    
    numericVals <- toPlot$value[!is.infinite(toPlot$value)]
    minNumeric <- min(numericVals) - 2.5
    maxNumeric <- max(numericVals) + 2.5
    maxAbs <- max(abs(c(minNumeric, maxNumeric)))
    minAbs <- min(abs(c(minNumeric, maxNumeric)))
    toPlot$value[is.infinite(toPlot$value) & toPlot$value < 0] <- -maxAbs
    toPlot$value[is.infinite(toPlot$value) & toPlot$value > 0] <- maxAbs
    
    p <- ggplot(toPlot, aes(y = gene_name, x = value)) +
        # Points
        geom_point(aes(colour = type, group = type), position = position_dodgev(height = 0.3)) +
        scale_colour_manual(values = c("#0073e6", "#e65800")) +
        coord_cartesian(xlim = c(-maxAbs - 0.15, maxAbs + 0.15)) +
        geom_vline(xintercept = 0, colour = "grey30", linetype = "dotted") +
        # Mean stripes
        annotate("rect", colour = "#8AC4FF", fill = "#8AC4FF", xmin = dnds_mis_mean$conf_lo, xmax = dnds_mis_mean$conf_up, ymin = -Inf, ymax = Inf) +
        #annotate("text", colour = "#0073e6", x = dnds_mis_mean$mean, y = inf, label = "genome-wide mean", hjust = 1, vjust = 1) +
        #annotate("text", colour = "#0073e6", x = dnds_non_mean$mean, y = inf, label = "Nonsense mean", hjust = 1, vjust = 1) +
        annotate("rect", colour = "#FFC39E", fill = "#FFC39E", xmin = dnds_non_mean$conf_lo, xmax = dnds_non_mean$conf_up, ymin = -Inf, ymax = Inf) +
        geom_vline(xintercept = dnds_mis_mean$mean, colour = "#0073e6", linetype = "dashed") +
        geom_vline(xintercept = dnds_non_mean$mean, colour = "#e65800", linetype = "dashed") +
        # Inf stripes
        annotate("text", label = "-Inf", x = -Inf, y = Inf, hjust = 0, vjust = 1) +
        annotate("text", label = "Inf", x = Inf, y = Inf, hjust = 1, vjust = 1) +
        annotate("rect", colour = "grey30", alpha = 0.3, xmin = -Inf, xmax = -maxAbs + 1, ymin = -Inf, ymax = Inf) +
        annotate("rect", colour = "grey30", alpha = 0.3, xmax = Inf, xmin = maxAbs - 1, ymin = -Inf, ymax = Inf) +
        xlab("dN/dS") +
        ylab("") +
        theme_bw() + 
        theme(legend.position = "top")
    
    return(p)
}
        
    
plot_gene_presence_heatmap <- function(mutations, expected_mutations, normalize_expected = F) {
    
    maxMuts <- 20
    
    heatmap_cols <- colorRampPalette(c("#cee0ff", "#428aff"))(maxMuts + 1)
    #heatmap_cols <- c("grey", heatmap_cols)

    gene_count_per_tissue <- mutations %>%
        group_by(tissue, gene_name) %>%
        summarise(n_muts = n()) %>%
        ungroup() %>%
        spread(tissue, n_muts) 
    
    
    # Get gene order
    genes_ordered <- get_gene_orders(mutations)
    
    # Get tissue order
    gene_per_tissue <- mutations %>%
        group_by(tissue) %>%
        summarise(n_genes = length(unique(gene_name))) %>%
        ungroup()
    
    tissues_ordered <- gene_per_tissue$tissue[order(-gene_per_tissue$n_genes)]
        
    
    # Converting to numeric matrix
    gene_count_per_tissue <- as.data.frame(gene_count_per_tissue)
    rownames(gene_count_per_tissue) <- gene_count_per_tissue$gene_name
    gene_count_per_tissue <- gene_count_per_tissue[,-1]
    gene_count_per_tissue <- as.matrix(gene_count_per_tissue)
    
    # Normalize by total number of mutations per tissue
    if(normalize_expected) {
        
        # Get total number of predicted mutations per tissue and use it to normalize all the genes in that tissue
        for(i in 1:ncol(gene_count_per_tissue)) {
            tissue <- colnames(gene_count_per_tissue)[i]
            gene_count_per_tissue[,i] <-  gene_count_per_tissue[,i] / expected_mutations[tissue, "observed"]
        }
        
        # Order tissues based on total mutaiton rate
        tissues_ordered <- colnames(gene_count_per_tissue)[order(-colSums(gene_count_per_tissue, na.rm = T))]
        gene_count_per_tissue <- log10(gene_count_per_tissue)
        gene_count_per_tissue[is.infinite(gene_count_per_tissue)] <- NA
    }
    
    # Removing NAs and big numbers
    #gene_count_per_tissue[is.na(gene_count_per_tissue)] <- 0
    gene_count_per_tissue[gene_count_per_tissue > maxMuts ] <- maxMuts
    gene_count_per_tissue <- gene_count_per_tissue[genes_ordered,tissues_ordered]
    
    heatmap.2(gene_count_per_tissue, Rowv = F, Colv = F, trace = "none", col = heatmap_cols, 
              density.info = "none", key.xlab = "Mutations", symbreaks = F,
              colsep = 1:ncol(gene_count_per_tissue), rowsep = 1:nrow (gene_count_per_tissue), na.color = "grey",
              sepwidth=c(0.001,0.001), sepcolor = "black")
    
    return(invisible(NULL))
}

get_gene_orders <- function(mutations) {
    
    # Get gene order
    gene_count <- mutations %>%
        group_by(gene_name) %>%
        summarise(mut_count = n()) %>%
        ungroup()
    
    genes_ordered <- gene_count$gene_name[order(-gene_count$mut_count, gene_count$gene_name)]
    
    return(genes_ordered)
}

make_oncokb <- function(mutations) {
    data.frame(Hugo_Symbol =  mutations$gene_name,
               Entrez_Gene_Id = mutations$gene_id,
               Center = ".", 
               NCBI_Build = ".", 
               Chromosome =  mutations$chr,
               Start_Position =  mutations$pos,
               End_Position = mutations$pos,
               Strand =  "+",
               Variant_Classification = mutations$impact, 
               Variant_Type = "SNP", 
               Reference_Allele =  mutations$ref, 
               Tumor_Seq_Allele1 = mutations$ref, 
               Tumor_Seq_Allele2 =  mutations$mut, 
               dbSNP_RS =  ".", 
               dbSNP_Val_Status =  ".",
               Tumor_Sample_Barcode =  mutations$sample,
               Matched_Norm_Sample_Barcode =  mutations$sample,
               Match_Norm_Seq_Allele1 =  mutations$ref,
               Match_Norm_Seq_Allele2 =  mutations$ref,
               Tumor_Validation_Allele1 =  mutations$ref,
               Tumor_Validation_Allele2 = mutations$mut,
               Match_Norm_Validation_Allele1 =  mutations$mut,
               Match_Norm_Validation_Allele2 =  mutations$mut,
               Verification_Status =  "validated",
               Validation_Status = ".",
               Mutation_Status = "somatic",
               Sequencing_Phase = ".",
               Sequence_Source =  "GTEx",
               Validation_Method =  ".",
               Score =  ".",
               BAM_File =  ".",
               Sequencer =  ".",
               Tumor_Sample_UUID = ".", 
               Matched_Norm_Sample_UUID =  ".",
               HGVSc =  paste0("c.", gsub("[ATCGN]+", "", mutations$ntchange), mutations$ref_cod, ">", mutations$mut_cod),
               HGVSp =  ".",
               HGVSp_Short =  paste0("p.", mutations$aachange),
               Transcript_ID =  ".",
               Exon_Number =  ".",
               t_depth = mutations$coverage,
               t_ref_count = mutations$coverage - mutations$n_support,
               t_alt_count =  mutations$n_support, 
               n_depth =  mutations$coverage,
               n_ref_count = mutations$coverage - mutations$n_support,
               n_alt_count =  mutations$n_support,
               Consequence =  mutations$impact,
               Protein_position =  gsub("[A-Z\\*]+", "", mutations$aachange)
               )
}

main()    

# plots the outputs of mutations_RNAexpression_association.R
# This focuses on individual genes:
#   - Most mutated
#   - Most tissue-spefific

library("gplots")
library("tidyr")
library("dplyr")
library("ggplot2")
library("metap")
source("../../R/ggthemes.R", chdir = T)
source("../../R/mutationGeneAnnotation.R", chdir = T)

PRUNE_PVALS <- 10 # -log10
PVAL_SIGNIF <- 0.01

main <- function(cmdArgs = commandArgs(T)) {
    
    out_prefix <- cmdArgs[1]
    expression_association_folder <- cmdArgs[2]
    
    #out_prefix <- "~/deleteme_out_expre"
    #expression_association_folder <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations"
    
    # Read tables
    expression_association <- read_files(expression_association_folder)
    
    
    # Make heatmaps
    for(mut_type in unique(expression_association$snps)) {

        # Select current mutation
        current_association <- expression_association [expression_association$snps == mut_type,]
        
        # Plot general numbers of significant genes across tissues
        plot_gene_tissue_numbers(current_association, paste0(out_prefix, "gene_plots_", mut_type))
        
    }
    
    writeLines("done", paste0(out_prefix, "plots_done_genes.txt"))
    
}
    

read_files <- function(x) {
    
    results <- list()
    tissue_folders <- list.dirs(x, recursive = F)
    for(i in tissue_folders) {
        result_files <- list.files(i)[grepl("ignificant.txt", list.files(i))]
        for(j in result_files) {
            flush.console()
            cat(i, "...", j, "\n")
            current <- read.table(gzfile(file.path(i,j)), sep = "\t", stringsAsFactors = F, header = T)
            current$tissue <- basename(i)
            results[[paste(i,j)]] <- current
        }
    }
    
    results <- do.call(rbind, results)
    rownames(results) <- 1:nrow(results)
    return(results)
}

# Takes gene association results of one mutation type and 
# plots the following:
#  - Histogram of number of tissues a gene is significant on
#  - Table with how many genes per bin ther are from the above histogram
#  - Histogram of number of tissues a gene is significant on 
#  - Histogram percentage of tissues (from total tested)a gene is significant on
plot_gene_tissue_numbers <- function(x, out_prefix) {
    
    # Counts how many significant associations (tissues) there are per gene
    signif_per_gene <- x %>%
        group_by(gene) %>%
        summarise(significant = sum(pvalue_bonf < PVAL_SIGNIF), non = sum(pvalue_bonf >= PVAL_SIGNIF), total = significant + non, p_fisher = ifelse(n() > 1, sumlog(pvalue_bonf)$p, pvalue_bonf),
                  significant_pos = sum(pvalue_bonf < PVAL_SIGNIF & statistic > 0), significant_neg = sum(pvalue_bonf < PVAL_SIGNIF & statistic < 0),
                  p_fisher_pos =  ifelse(sum(statistic > 0) == 0, 1, ifelse(sum(statistic > 0) > 1,  sumlog(pvalue_bonf[statistic > 0])$p, pvalue_bonf[statistic > 0])),
                  p_fisher_neg =  ifelse(sum(statistic < 0) == 0, 1, ifelse(sum(statistic < 0) > 1,  sumlog(pvalue_bonf[statistic < 0])$p, pvalue_bonf[statistic < 0]))
                  ) %>%
        ungroup() %>%
        mutate(ratio = significant/(significant + non), ratio_pos = significant_pos/(significant + non), ratio_neg= significant_neg/(significant + non))
    
    # Histograms of number of significant tissues per gene and number of total tissues per gene tests
    p_hist_genes_tissue <- plot_hist_percentages(signif_per_gene, "significant") + xlab("Tissues a gene is significant on") + theme_grid_y()
    p_hist_genes_tested <- plot_hist_percentages(signif_per_gene, "total", totalLabel = F, round_percent_n = 0) + xlab("Tissues a gene was tested on") + theme_grid_y()
    
    # Histogram of percentage of tissues where a gene is significant from all tests only for 
    # genes that were teste for a high number of tissues (above median numbe of tests)
    median_n_tests <- median(signif_per_gene$total)
    signif_per_gene_high_tests <- signif_per_gene[signif_per_gene$total >= median_n_tests,]
    p_hist_percentage_pos_highly_tested <- plot_hist_percentages(signif_per_gene_high_tests, "ratio", n_bins = 15) + 
                                            xlab("Percentage of tissues a gene is significant on")+
                                            theme_grid_y() +
                                            theme(axis.text.x = element_text(angle = 30, hjust = 1))
    
    # Plot heatmap of most shared genes
    signif_per_gene_high_tests <- signif_per_gene_high_tests[order(-signif_per_gene_high_tests$ratio, signif_per_gene_high_tests$p_fisher),]
    most_shared_genes <- head(signif_per_gene_high_tests$gene, 40)
    most_shared <- x[x$gene %in% most_shared_genes, ]
    most_shared$signed_pvalue <- -log10(most_shared$pvalue_bonf)
    most_shared$signed_pvalue[most_shared$statistic < 0] <- -most_shared$signed_pvalue[most_shared$statistic < 0] 
    
    plot_pvalue_heatmap(most_shared, gene, tissue, signed_pvalue, out_prefix)
    
    # Write genes table ordered by in how many tissues it is significant
    writeLines(signif_per_gene_high_tests$gene, paste0(out_prefix, "genes_ordered_by_number_of_tissues_is_significant.txt"))
    
    # Write genes table ordered by in how many tissues they are possitively associated
    signif_per_gene_high_tests <- signif_per_gene_high_tests[order(-signif_per_gene_high_tests$ratio_pos, signif_per_gene_high_tests$p_fisher_pos),]
    writeLines(signif_per_gene_high_tests$gene, paste0(out_prefix, "genes_ordered_by_number_of_tissues_is_positively_significant.txt"))
    
    # Write genes table ordered by in how many tissues they are possitively associated
    signif_per_gene_high_tests <- signif_per_gene_high_tests[order(-signif_per_gene_high_tests$ratio_neg, signif_per_gene_high_tests$p_fisher_neg),]
    writeLines(signif_per_gene_high_tests$gene, paste0(out_prefix, "genes_ordered_by_number_of_tissues_is_negatively_significant.txt"))
    
    # Save other plots
    ggsave(paste0(out_prefix, "_histogram_signif_genes_per_tissue.pdf"), p_hist_genes_tissue, height = 4, width = 9)
    ggsave(paste0(out_prefix, "_histogram_tests_per_tissue.pdf"), p_hist_genes_tested, height = 4, width = 9)
    ggsave(paste0(out_prefix, "_histogram_percentage_tissue_positive_highly_tested.pdf"), p_hist_percentage_pos_highly_tested, height = 4, width = 9)
    
}

plot_hist_percentages <- function(x, column, percentageLabel = T, totalLabel = T, round_percent_n = 2, n_bins = NULL) {
    
    if(is.null(n_bins)) {
        to_plot <- as.data.frame(table(x[,column]))
    } else {
        bins <- seq(min(x[,column]), max(x[,column]), length.out = n_bins)
        binsLabels <- signif(bins*100,2)
        binsLabels <- paste0(binsLabels[-n_bins], " - ", binsLabels[-1])
        countBins <- cut(x[,column, drop = T], breaks = bins, labels = binsLabels, include.lowest = T)
        to_plot <- as.data.frame(table(countBins))
    }
    
    colnames(to_plot) <- c("Var1", "Freq")
    to_plot$cum_total <- cumsum(to_plot$Freq)
    to_plot$cum_percentage <- to_plot$cum_total / sum(to_plot$Freq) * 100
    to_plot$label <- ""
    
    if(percentageLabel)
        to_plot$label <- paste0(to_plot$label, round(to_plot$cum_percentage, round_percent_n), "%")
    if(totalLabel)
        to_plot$label <- paste0(to_plot$label, "\nn = ", to_plot$Freq)
    
    p <- ggplot(to_plot, aes(x = Var1, y = Freq)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = label), vjust = 0)
        
   return(p)
}


plot_pvalue_heatmap <- function(x, variable1, variable2, value, out_prefix) {
    
    variable1<- enquo(variable1)
    variable2<- enquo(variable2)
    value<- enquo(value)
    
    # Do heatmap
    # Spread into matrix of pvalues with dimensions gene by tissues
    current_association_matrix <- spread_tissues(x, variable1, variable2, value)
    current_association_matrix[current_association_matrix > PRUNE_PVALS] <- PRUNE_PVALS
    current_association_matrix[current_association_matrix < -PRUNE_PVALS] <- -PRUNE_PVALS
    # Order rows by number of significant interactions
    temp_mat <- current_association_matrix
    temp_mat[is.na(temp_mat)] <- 0
    
    current_association_matrix <- current_association_matrix[order(-rowSums(abs(temp_mat) > -log10(PVAL_SIGNIF))),]
    current_association_matrix <- current_association_matrix[,order(-colSums(abs(temp_mat) > -log10(PVAL_SIGNIF)))]
    
    rownames(current_association_matrix) <- paste0(rownames(current_association_matrix), " - ", queryRefCDS(rownames(current_association_matrix), reference = "gene_id", query = "gene_name"))
    
    heatmap.pars <- list(trace = "none", density.info = "none", Rowv = F, dendrogram = "none", Colv = F,
                         symkey = T, col = get_colors_significant(PRUNE_PVALS, -log10(PVAL_SIGNIF))
                         )
    
    pdf(paste0(out_prefix, "_heatmap_associations_", ".pdf"), height = 6, width = 9)
    do.call(heatmap.2, c(list(x=current_association_matrix),
                         heatmap.pars))
    dev.off()
}

# From a data frame containing the cols alias, signed_pvalue, tissue
# spreads signed_pvalues into tissue columns
spread_tissues <- function(x, variable1, variable2, value) {
    
    x <- x %>%
      dplyr::select(!!variable1, !!value, !!variable2) %>%
      spread(!!variable2, !!value)
    
    x <- as.data.frame(x)
    rownames(x) <- x[,1]
    x <- x[,-1]
    
    return(as.matrix(x))
}

# Create color palette for only significant based on log10 pvalue
get_colors_significant <- function(maxPval = 10, significant = 2) {
    
    blues <- colorRampPalette(c("#4949FF", "#2F2F50"))(maxPval - significant + 1)
    yellows <- colorRampPalette(c("#50502F", "#FFFF4A"))(maxPval - significant + 1)
    blacks <- rep("grey20", significant * 2 -1)
    
    return(c(blues,blacks,yellows))
}

main()

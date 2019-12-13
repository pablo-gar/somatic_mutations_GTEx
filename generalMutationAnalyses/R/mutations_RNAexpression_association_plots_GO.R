# plots the outputs of mutations_RNAexpression_association.R
# This focuses on GO categories
#   - Most mutated
#   - Most tissue-spefific

library("gplots")
library("tidyr")
library("dplyr")
library("ggplot2")
source("../../R/ggthemes.R", chdir = T)
source("../../R/mutationGeneAnnotation.R", chdir = T)

PRUNE_PVALS <- 10 # -log10
PVAL_SIGNIF <- 0.01
GO_TOP_N <- 10

GO_TERMS <- c("vir", "immu", "repair|damage", "local|transp", "neur|nervo")


main <- function(cmdArgs = commandArgs(T)) {
    
    out_prefix <- cmdArgs[1]
    expression_association_folder <- cmdArgs[2]
    
    #out_prefix <- "out_"
    #expression_association_folder <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations"
    
    # Read tables
    expression_association <- read_files(expression_association_folder, "ignificant.txt")
    
    # Read GO results
    GO_results <- read_files(expression_association_folder, "GO.txt$")
    
    # Plot top GO results per tissue per mutation type
    for(mut_type in unique(GO_results$mut_type)) {
        
        # Select current mutation
        current_GO <- GO_results [GO_results$mut_type == mut_type,]
        
        # Select significant terms
        current_GO_signif <- current_GO[current_GO$p_bonf < PVAL_SIGNIF,]
        
        # Plot top 10 categories per tissue
        plot_top_GO(current_GO_signif, n_GO = GO_TOP_N, out_prefix = paste0(out_prefix, mut_type))
        
        # Plot virus related GOs
        forPlot <- current_GO_signif[grep("vir", current_GO_signif$Term, perl = T),]
        if(nrow(forPlot) > 0) 
            plot_top_GO(forPlot, n_GO = 1e5, out_prefix = paste0(out_prefix, "virus_", mut_type))
        
    }
    
    # Plot specific GO terms
    for(go_term in GO_TERMS) {
        
        # Select current mutation
        current_GO <- GO_results[grep(go_term, GO_results$Term, perl = T),]
        
        # Select significant terms
        current_GO_signif <- current_GO[current_GO$p_bonf < PVAL_SIGNIF,]
        
        # plot
        if(nrow(current_GO_signif) > 0) {
            current_GO_signif$tissue <- paste0(current_GO_signif$mut_type, "-", current_GO_signif$tissue)
            current_GO_signif <- current_GO_signif[current_GO_signif$mut_type != "all",]
            plot_top_GO(current_GO_signif, n_GO = 15, out_prefix = paste0(out_prefix, go_term ))
        }
        
    }
    
    writeLines("done", paste0(out_prefix, "plots_done_GO.txt"))
    
}
    

read_files <- function(x, suffix) {
    
    results <- list()
    tissue_folders <- list.dirs(x, recursive = F)
    for(i in tissue_folders) {
        result_files <- list.files(i)[grepl(suffix, list.files(i), perl = T)]
        for(j in result_files) {
            flush.console()
            cat(i, "...", j, "\n")
            current <- read.delim(file.path(i,j), sep = "\t", stringsAsFactors = F, header = T)
            current$tissue <- basename(i)
            results[[paste(i,j)]] <- current
        }
    }
    
    results <- do.call(rbind, results)
    rownames(results) <- 1:nrow(results)
    return(results)
}


plot_top_GO <- function(x, n_GO, out_prefix) {
    
    nCol <- 3
    width_panel <- 7
    height_panel <- 0.4
    width_plot <- width_panel * nCol
    
    # Get ranks
    x <- x %>%
        group_by(tissue) %>%
        mutate(rank = rank(p_bonf, ties.method = "first"), log10_pval = -log10(p_bonf)) %>%
        ungroup()
    
    # only top 10 terms
    x <- x[x$rank <= n_GO,]
    
    # Order terms
    x <- x[order(x$tissue, x$rank, x$Term),]
    tissues <- setNames(1:length(unique(x$tissue)), unique(x$tissue))
    x$Term <- paste0(x$Term, "-", tissues[x$tissue])
    x$Term <- factor(x$Term, levels = rev(unique(x$Term)), ordered = T)
    
    # max number of rows per panel
    rows_per_panel <- max(x$rank)
    
    # Plot
    p <- ggplot(x, aes(y = Term, x = log10_pval)) +
        geom_segment(aes(yend = Term), colour = "grey50", xend = 0) +
        geom_point(size = 2) +
        facet_wrap(~tissue, scales = "free_y", ncol = nCol) +
        theme_grid_x()
    
    ggsave(paste0(out_prefix, "_topGOresults_per_tissue.pdf"), p, width = width_plot, height = ceiling(length(tissues)/ nCol)  * height_panel * rows_per_panel, limitsize = FALSE)
    

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
    
    pdf(paste0(out_prefix, "_heatmap_associations_", mut_type, ".pdf"), height = 6, width = 9)
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

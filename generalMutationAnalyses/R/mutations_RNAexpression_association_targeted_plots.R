# plots the outputs of mutations_RNAexpression_association_targeted.R
# it divides genes in type of repair and plots the associations pvalues in heatmaps
library("gplots")
library("tidyr")
library("dplyr")
library("metap")
library("ggplot2")
source("../../R/ggthemes.R", chdir = T)
source("../../R/mutationGeneAnnotation.R", chdir = T)
source("../../R/FDR.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    PRUNE_PVALS <- 10 # -log10
    PVAL_SIGNIF <- 0.01
    
    pathway_gene_file <- cmdArgs[1]
    out_prefix <- cmdArgs[2]
    workingDir <- cmdArgs[3]
    
    #pathway_gene_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/repair_genes.txt"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations_targeted_plots/n6_0.0_0.7_"
    #workingDir <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations_targeted"
    
    # Gathering the three different types of files
    expression_association_files <- list.files(workingDir, recursive = T, full.names = T, pattern = ".*\\-results.txt")
    gene_group_pvalues_files <- list.files(workingDir, recursive = T, full.names = T, pattern = ".*_group_enrichment.txt")
    permuted_tables_files <- list.files(workingDir, recursive = T, full.names = T, pattern = ".*_permuted_tables.Rdata")
    
    # Read tables
    expression_association <- read_files(expression_association_files)
    rownames(expression_association) <- paste0(expression_association$alias, ".", expression_association$tissue, ".", expression_association$snps)

    permuted_tables <- read_permuted_tables(permuted_tables_files)
    permuted_tables$real_p <- expression_association[ paste0(permuted_tables$gene, ".", permuted_tables$tissue, ".", permuted_tables$mut_type), "pvalue"]
    nPerm <- max(permuted_tables$permutation)
    
    pathway_gene_table <- read.table(pathway_gene_file, sep = "\t", stringsAsFactors = F, row.names = 1, header = T)
    gene_group_pvalues <- read_files(gene_group_pvalues_files)
    
    
    # Convert p-values to signed log10 pvalues
    expression_association <- read_files(expression_association_files)
    expression_association$signed_pvalue <- -log10(expression_association$pvalue_bonf) * ifelse(expression_association$statistic > 0, 1, -1)
    
    # Adding pathway name per gene
    pathway_gene_table <- pathway_gene_table[rownames(pathway_gene_table) %in% expression_association$alias,, drop =F]
    expression_association$pathway <- pathway_gene_table[expression_association$alias, "pathway"]
    
    # Make heatmaps
    for(mut_type in unique(expression_association$snps)) {
        
        # Select current mutation
        current_association <- expression_association[expression_association$snps == mut_type, ]
        
        ##### bonf pvalues of the linear regression 
        
        # Do bar plot of number of positive and negative interactions
        significant_per_gene <- get_counts_per_gene(current_association, pval_cutoff = PVAL_SIGNIF, pathway_gene_table = pathway_gene_table)
        this_out_prefix <- paste0(out_prefix, "association_barcount_", mut_type, ".pdf")
        do_barplot(significant_per_gene, this_out_prefix)
        
        # Do heatmap
        this_out_prefix <- paste0(out_prefix, "heatmap_associations_", mut_type, ".pdf")
        tissueOrder <- do_heatmap(current_association, pathway_gene_table, PRUNE_PVALS, PVAL_SIGNIF, this_out_prefix) 
        
        
        # Plot gene enrichment
        gene_enrichment_p <- get_gene_enrichment(permuted_tables, mut_type)
        this_out_prefix <- paste0(out_prefix, "heatmap_associations_gene_", mut_type, ".pdf")
        pdf(this_out_prefix, width = 4)
        do_heatmap_pvalues(gene_enrichment_p[rownames(pathway_gene_table), c(1,1)])
        dev.off()
        
        # Plot pathways enrichment
        pathway_enrichment_p <- get_pathway_erichment(gene_group_pvalues, mut_type)
        this_out_prefix <- paste0(out_prefix, "heatmap_associations_groups_", mut_type, ".pdf")
        pdf(this_out_prefix, height = 6)
        do_heatmap_pvalues(pathway_enrichment_p[unique(pathway_gene_table$pathway), tissueOrder])
        dev.off()
        
        
        # Get FDR enrichment per gene
        genes <- unique(expression_association$alias)
        gene_fdrs <- get_gene_fdrs(expression_association, mut_type, genes, permuted_tables)
        this_out_prefix <- paste0(out_prefix, "heatmap_associations_genes_fdr_", mut_type, ".pdf")
        pdf(this_out_prefix, height = 6)
        do_heatmap_fdrs(gene_fdrs[rownames(pathway_gene_table), c(1,1)])
        dev.off()
        write.table(gene_fdrs, paste0(this_out_prefix, ".txt"), sep = "\t", quote = F, col.names = T, row.names = T)
        
        # Get FDR enrichment per pathway
        groups <- by(rownames(pathway_gene_table), pathway_gene_table$pathway, function(x) as.character(x))
        pathway_fdrs <- get_pathway_fdrs(expression_association, mut_type, groups, permuted_tables)
        this_out_prefix <- paste0(out_prefix, "heatmap_associations_groups_fdr_", mut_type, ".pdf")
        pdf(this_out_prefix, height = 6)
        do_heatmap_fdrs(pathway_fdrs[unique(pathway_gene_table$pathway), tissueOrder])
        dev.off()
        
        this_out_prefix <- paste0(out_prefix, "heatmap_associations_groups_fdr_clustered_", mut_type, ".pdf")
        pdf(this_out_prefix, height = 6)
        do_heatmap_fdrs(pathway_fdrs[unique(pathway_gene_table$pathway), tissueOrder], Colv = T)
        dev.off()
        
        # write table fdr pathways
        pathway_fdrs <- as.data.frame(pathway_fdrs)
        pathway_fdrs$pathway <- rownames(pathway_fdrs)
        pathway_fdrs <- gather(pathway_fdrs, "tissue", "fdr", -pathway)
        write.table(pathway_fdrs, paste0(this_out_prefix, ".txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        
        
        
        ######################################### 
        ##### permuted pvalues -- Gene level permutations
        #current_association <- current_association %>%
        #    group_by(tissue) %>%
        #    mutate(pvalue_perm_bonf = p.adjust(pvalue_perm, method = "BH") - 0.00001, signed_pvalue_perm = -log10(pvalue_perm_bonf) * ifelse(statistic > 0, 1, -1)) %>%
        #    ungroup() %>%
        #    as.data.frame() 
        #        
        #current_association$signed_pvalue <- current_association$signed_pvalue_perm
        #current_association$pvalue_bonf <- current_association$pvalue_perm_bonf
        #
        ## Do bar plot of number of positive and negative interactions
        #significant_per_gene <- get_counts_per_gene(current_association, pval_cutoff = 0.05, pathway_gene_table = pathway_gene_table)
        #this_out_prefix <- paste0(out_prefix, "association_barcount_permPval", mut_type, ".pdf")
        #do_barplot(significant_per_gene, this_out_prefix)
        #
        ## Do heatmap
        #this_out_prefix <- paste0(out_prefix, "heatmap_associations_permPval", mut_type, ".pdf")
        #do_heatmap(current_association, pathway_gene_table, 3, 0.05, this_out_prefix) 
    }

    writeLines("done", paste0(out_prefix, "done_targeted_plots.txt"))
    
}


read_files <- function(x) {
    
    results <- list()
    for(i in x) {
        current <- read.table(i, sep = "\t", stringsAsFactors = F, header = T)
        current$tissue <- basename(dirname(i))
        results[[i]] <- current
    }
    
    results <- do.call(rbind, results)
    rownames(results) <- 1:nrow(results)
    return(results)
}

read_permuted_tables <- function(x) {
    
    results <- list()
    counter <- 1
    for(i in x) {
        tissue <- basename(dirname(i))
        load(i)
        for(mut_type in names(x_perm_all)) {
            current <- as.data.frame(x_perm_all[[mut_type]])
            current$permutation <- 1:nrow(current)
            current <- gather(current, "gene", "pvalue", -permutation)
            current$gene <-  queryRefCDS(current$gene, reference = "gene_id", query = "gene_name")
            current$mut_type <- mut_type
            current$tissue <- tissue
            results[[counter]] <- current
            counter <- counter + 1
        }
    }
    
    results <- do.call(rbind, results)
    
    return(results)
    
}

# From a data frame containing the cols alias, signed_pvalue, tissue
# spreads signed_pvalues into tissue columns
spread_tissues <- function(x, column) {
    
    column_en <- enquo(column)
    
    x <- x %>%
      dplyr::select(alias, !!column_en, tissue) %>%
      spread(tissue, !!column_en)
    
    x <- as.data.frame(x)
    rownames(x) <- x[,1]
    x <- x[,-1]
    
    return(as.matrix(x))
}
 
# From a data frame containing the cols stastic, pval_bonf, alias
# counts the number of positive and negative signifcant hits
get_counts_per_gene <- function(x, pval_cutoff, pathway_gene_table) {
    x <- x %>%
      group_by(alias) %>%
      filter(pvalue_bonf < pval_cutoff) %>%
      summarise(Positive = sum(statistic > 0), Negative = sum(statistic < 0)) %>%
      ungroup() %>%
      gather("Association_direction", "count", Positive, Negative)
  
  x$pathway <- factor(pathway_gene_table[x$alias, "pathway"], levels = unique(pathway_gene_table$pathway), ordered = T)
  x$alias <- factor(x$alias, levels = rev(rownames(pathway_gene_table)), ordered = T)
  return(x)
}

# Create color palette for only significant based on log10 pvalue
get_colors_significant <- function(maxPval = 10, significant = 2) {
    
    added <- (significant %% 1 > 0) * 2
    
    blues <- colorRampPalette(c("#4949FF", "#2F2F50"))(maxPval - floor(significant))
    yellows <- colorRampPalette(c("#50502F", "#FFFF4A"))(maxPval - floor(significant))
    blacks <- rep("grey20", floor(significant) * 2 -1 + added)
    
    return(c(blues,blacks,yellows))
}

get_breaks <- function(maxPval = 10, significant = 2) {
    
    pos_breaks <- unique(sort(c(1:maxPval, significant)))
    
    return(c(rev(-pos_breaks), pos_breaks))
}

do_barplot <- function(significant_per_gene, this_out_prefix) {
    
    p_bars <- ggplot(significant_per_gene, aes(x = alias, y = count)) +
        geom_bar(aes(fill = Association_direction), stat = "identity") +
        ylab("No. of significant associations") +
        scale_fill_manual(values = c("#c65353", "#6699ff")) +
        #facet_wrap(.~pathway, scales = "free_y", ncol = 1) +
        #facet_grid(~pathway, scales = "free_y") +
        theme_grid_x() +
        coord_flip() +
        theme(legend.position = "top")
    
    ggsave(this_out_prefix, p_bars, width = 4)
    
}

do_heatmap <- function(current_association, pathway_gene_table, prune_pvals, pval_signif, this_out_prefix) {
    
    # Spread into matrix of pvalues with dimensions gene by tissues
    current_association_matrix <- spread_tissues(current_association, signed_pvalue)
    current_association_matrix <- current_association_matrix[rownames(pathway_gene_table),]
    current_association_matrix[current_association_matrix > prune_pvals] <- prune_pvals
    current_association_matrix[current_association_matrix < -prune_pvals] <- -prune_pvals
    
    
    
    #Gets fdr in matrix form
    current_association_matrix_fdr <- spread_tissues(current_association, fdr)
    current_association_matrix_fdr <- current_association_matrix_fdr[rownames(pathway_gene_table),]
    
    # Get tissue orders
    Colv <- hclust(dist(t(current_association_matrix)))
    tissues <- Colv$labels[Colv$order]
    Colv <- as.dendrogram(Colv)
    
    # Re order base on hclust
    current_association_matrix_fdr <- current_association_matrix_fdr[,tissues]
    current_association_matrix <- current_association_matrix[,tissues]
    
    
    # Get FDR labels
    fdr_labels <- current_association_matrix_fdr
    fdr_labels[,] <-  ""
    fdr_labels[current_association_matrix_fdr <= 0.2] <- "*"
    fdr_labels[current_association_matrix_fdr <= 0.1] <- "**"
    fdr_labels[current_association_matrix_fdr <= 0.05] <- "***"
    
    
    
    heatmap.pars <- list(trace = "none", density.info = "none", Rowv = F, Colv = F, dendrogram = "none", cellnote = fdr_labels,
                         symkey = T, col = get_colors_significant(prune_pvals, -log10(pval_signif)), breaks = get_breaks(prune_pvals, -log10(pval_signif))
                         )
    
    #heatmap.pars <- list(trace = "none", density.info = "none", Rowv = F, Colv = F, dendrogram = "none", cellnote = fdr_labels,
    #                     symkey = T, col =  colorRampPalette(c("#4949FF", "#FFF4AA"))(10)
    #                     )
    
    pdf(this_out_prefix, height = 5)
    do.call(heatmap.2, c(list(x=current_association_matrix),heatmap.pars))
    dev.off()
    
    
    return(tissues)
    
}

do_heatmap_pvalues <- function(x) {
    
    y <- x
    
    # Convert pvalues
    x[ y > 0.05 ] <- 0
    x[ y > 0.01 & y <= 0.05 ] <- 1
    x[ y > 0.001 & y <= 0.01 ] <- 2
    x[ y > 0.0001 & y <= 0.001 ] <- 3
    x[ y <= 0.0001 ] <- 4
    
    
    
    heatmap.pars <- list(x = x, trace = "none", density.info = "none", Rowv = F, Colv = F, dendrogram = "none",
                         col = colorRampPalette(c("white", "red"))(max(x) + 1),
                         colsep = 1:ncol(x),
                         rowsep = 1:nrow(x),
                         sepcolor = "grey30"
                         )
    
   # pdf(this_out_prefix, height = 5)
   do.call(heatmap.2, heatmap.pars)
   # dev.off()
    
}

do_heatmap_fdrs <- function(x, Colv = F, Rowv = F) {
    
    y <- x
    
    # Convert pvalues
    x[ y > 0.2 ] <- 0
    x[ y > 0.15 & y <= 0.20 ] <- 1
    x[ y > 0.1 & y <= 0.15 ] <- 2
    x[ y > 0.05 & y <= 0.1 ] <- 3
    x[ y > 0.01 & y <= 0.05 ] <- 4
    x[ y > 0.001 & y <= 0.01 ] <- 5
    x[ y > 0.0001 & y <= 0.001 ] <- 6
    x[ y <= 0.0001 ] <- 7
    
    
    
    heatmap.pars <- list(x = x, trace = "none", density.info = "none", Rowv = Rowv, Colv = Colv, dendrogram = "none",
                         col = colorRampPalette(c("white", "red"))(max(x) + 1),
                         colsep = 1:ncol(x),
                         rowsep = 1:nrow(x),
                         sepcolor = "grey30"
                         )
    
   # pdf(this_out_prefix, height = 5)
   do.call(heatmap.2, heatmap.pars)
   # dev.off()
    
}


get_gene_enrichment <- function(permuted_tables, mut_type) {
    
    current_perms <- permuted_tables[permuted_tables$mut_type == mut_type, ]
    
    gene_pvalue <- current_perms %>%
        group_by(gene, permutation) %>%
        summarise(fisher_permuted = sumlog(pvalue)$p, fisher_real = sumlog(real_p)$p) %>%
        ungroup() %>%
        group_by(gene) %>%
        summarise(enrichment_p = sum(fisher_permuted <= fisher_real)/n()) %>%
        ungroup() %>%
        as.data.frame()
    
    rownames(gene_pvalue) <- gene_pvalue[,1]
    gene_pvalue <- as.matrix(gene_pvalue[,-1, drop = F])
    
    return(gene_pvalue)
        
}

get_pathway_erichment <- function(gene_group_pvalues, mut_type) {
    
    gene_group_pvalues <- gene_group_pvalues[gene_group_pvalues$mut_type == mut_type, ]

    pathway_enrichment_p <- gene_group_pvalues %>%
           dplyr::select(-mut_type) %>%
           spread(tissue, pvalue) %>%
           as.data.frame()

    rownames(pathway_enrichment_p) <-  pathway_enrichment_p[,1]
    pathway_enrichment_p  <- as.matrix(pathway_enrichment_p[,-1])
    
}

#' @param groups is a named list
get_pathway_fdrs <- function(expression_association, mut_type, groups, permuted_tables) {
    
    tissues <-  unique(expression_association$tissue)
    results <- matrix(1, nrow = length(groups), ncol = length(tissues), dimnames = list(names(groups), tissues))
    permuted_tables <- permuted_tables[permuted_tables$mut_type == mut_type,]
    expression_association <- expression_association[expression_association$snps == mut_type,]
    
    for(tissue in tissues) {
        for(i in seq_along(groups)) {
            genes <- groups[[i]]
            group <- names(groups)[i]
            real_p <- expression_association[expression_association$alias %in% genes & expression_association$tissue == tissue, "pvalue"]
            random_p <- permuted_tables[ permuted_tables$gene %in% genes & permuted_tables$tissue == tissue, "pvalue"]
            min_fdr <- min(fdr_group(real_p, random_p))
            results[group, tissue] <- min_fdr
        }
    }
    
    return(results)

    
}

#' @param groups is a named list
get_gene_fdrs <- function(expression_association, mut_type, genes, permuted_tables) {
    
    results <- matrix(1, nrow = length(genes), ncol = 1, dimnames = list(genes, "all"))
    permuted_tables <- permuted_tables[permuted_tables$mut_type == mut_type,]
    expression_association <- expression_association[expression_association$snps == mut_type,]
    
    for(i in seq_along(genes)) {
        gene <- genes[i]
        real_p <- expression_association[expression_association$alias %in% gene, "pvalue"]
        random_p <- permuted_tables[ permuted_tables$gene %in% gene, "pvalue"]
        min_fdr <- min(fdr_group(real_p, random_p))
        results[gene, 1] <- min_fdr
    }
    
    return(results)

    
}


#' Performs the fdr values of enrichment at each pvalue cut off, with the pvalue
#' cutoffs being the pvalues of the members of the group. It uses a permuted
#' table that contains the same number of columns as pvalues and the number of
#' rows is the permutations
#' @param
fdr_group <- function(pvals, perm_matrix){
    
    observed_percentage <- rep(1, length(pvals)) 
    expected_percentage <- rep(1, length(pvals))
    
    for(i in seq_along(observed_percentage)) {
        expected_percentage[i] <- sum(perm_matrix <= pvals[i]) / length(perm_matrix)
        observed_percentage[i] <- sum(pvals <= pvals[i]) / length(pvals)
    }
    
    fdrs <- expected_percentage / observed_percentage
    fdrs[fdrs > 1] <- 1
    
    return(fdrs)
}
    
    

main()

# Gets GO enrichments of genes that show the highest strand bias across
# all mutation types

source("../../R/GO_analysis.R")
source("../../R/mutationGeneAnnotation.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    FDR_CUTOFF_STRAND_FC <- 0.05
    N_GO_CATEGORIES <- 50 # Categories to show
    TEST <- "fisher"
    
    out_prefix <- cmdArgs[1]
    fc_table_file <- cmdArgs[2]
    
    #out_prefix <- "out_"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Lung-n6_0.0_0.7.per_gene_fc_strand.txt"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Breast_Mammary_Tissue-n6_0.0_0.7.per_gene_fc_strand.txt"
    
    fc_table <- read.table(fc_table_file, sep = "\t", stringsAsFactors = F, header = T)
    fc_table$gene_id <- queryRefCDS(x = fc_table$gene, "gene_name", "gene_id")
    
    # Converting to log space
    fc_table$signed_log_fdr <- -log2(fc_table$fdr)
    fc_table$signed_log_fdr[fc_table$fc <= 1] <- -fc_table$signed_log_fdr[fc_table$fc <= 1]
    fc_table <- fc_table[!is.na(fc_table$signed_log_fdr),]
    
    for(mut_type in unique(fc_table$mutation)) {
        
        current <- fc_table[ fc_table$mutation == mut_type, ]
        
        
        # Perfomr GO analysis
        
        this_prefix <- paste0(out_prefix, "mut_type_", mut_type)
        perform_GO(current, this_prefix = this_prefix, fdr_fc = FDR_CUTOFF_STRAND_FC, n_GO = N_GO_CATEGORIES, statistic = TEST)
        
        
        #Saving gene tables
        write.table(current[current$fdr < 0.05 & current$fc >1,], paste0(this_prefix, "_high_fc_genes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        write.table(current[current$fdr < 0.05 & current$fc <= 1,], paste0(this_prefix, "_bottom_fc_genes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        
        # Save gene tables for GOrilla
        write.table(current[order(-current$fc),], paste0(this_prefix, "_high_fc_genes_GORILLA.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        write.table(current[order(current$fc),], paste0(this_prefix, "_bottom_fc_genes_GORILLA.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        
    }
    
    writeLines("Done", paste0(out_prefix, "allDone.txt"))
    
}
        
perform_GO <- function(current, this_prefix, fdr_fc, n_GO, statistic) {
    
    
    # Getting the top/bottom FC genes
    genes_test <- setNames(current$signed_log_fdr, current$gene_id)
    genes_test <- genes_test[!is.na(genes_test)]
    
    # GO enrichment
    if(statistic == "ks") {
        top_GO <- performGO(-genes_test, gene_selection_function = function(p) p, statistic = statistic)
        bottom_GO <- performGO(genes_test, gene_selection_function = function(p) p, statistic = statistic)
    } else {
        top_GO <- performGO(genes_test, gene_selection_function = function(p) p >= -log2(fdr_fc), statistic = statistic)
        bottom_GO <- performGO(genes_test, gene_selection_function = function(p) p <= log2(fdr_fc), statistic = statistic)
    }
    
    # Get GO results
    GO_tables_top <- get_GO_results(top_GO, test = statistic, fdr = F, bonf_signif = F, n = n_GO)
    GO_tables_bottom <- get_GO_results(bottom_GO, test = statistic, fdr = F, bonf_signif = F, n = n_GO)
    
    # Get genes associated with GO
    GO_genes_top <- get_genes_GO(genes = names(genes_test), x = GO_tables_top, GO_results = top_GO, n = n_GO)
    GO_genes_bottom <- get_genes_GO(genes = names(genes_test), x = GO_tables_bottom, GO_results = bottom_GO, n = n_GO)
    
    GO_genes_top$gene_name <- queryRefCDS(x = GO_genes_top$gene, "gene_id", "gene_name")
    GO_genes_bottom$gene_name <- queryRefCDS(x = GO_genes_bottom$gene, "gene_id", "gene_name")

    
    # Saving results
    write.table(GO_tables_top, paste0(this_prefix, "_high_fc_GO.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
    write.table(GO_tables_bottom, paste0(this_prefix, "_bottom_fc_GO.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
    
    write.table(GO_genes_top, paste0(this_prefix, "_high_fc_GO_genes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
    write.table(GO_genes_top, paste0(this_prefix, "_bottom_fc_GO_genes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
}

main()

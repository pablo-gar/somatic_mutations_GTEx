# takes the gene-level strand fc mutations from  plot_strandDifferences.R and
# calculates the correlations of fc

library("dplyr")
library("tidyr")
source("../../R/plots.R")
source("../../R/ggthemes.R")
#source("../../R/GO_analysis.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    file_1 <- cmdArgs[1]
    file_2 <- cmdArgs[2]
    tissue_1 <- cmdArgs[3]
    tissue_2 <- cmdArgs[4]
    out_prefix <- cmdArgs[5]
    
    #DEBUG
    #file_1 <- "/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences/Whole_Blood-n6_0.0_0.7.per_gene_fc_strand.txt"
    #file_2 <- "/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences_EXO/Whole_Blood_EXO-n6_0.0_0.7.per_gene_fc_strand.txt"
    #tissue_1 <- "RNA_seq"
    #tissue_2 <- "exome_seq"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/selection/validation_strand_fc_perGene_blood/"
    
    # Reading files
    fc_1 <- read.table(file_1, sep = "\t", stringsAsFactors = F, header = T)
    fc_2 <- read.table(file_2, sep = "\t", stringsAsFactors = F, header = T)
    fc_1$tissue <- tissue_1
    fc_2$tissue <- tissue_2
    
    # merging files
    merged_fc <- merge_tables(fc_1, fc_2)
    merged_fc <- merged_fc[ !is.na(merged_fc[,tissue_1]) & !is.na(merged_fc[,tissue_2]), ]
    
    # Perform scatter for all mutation types
    p <- scatter(merged_fc, x = tissue_1, y = tissue_2, facet_y = "mutation", ncol = 3, regression = T, method_cor = "sp", pSize = 1.5) +
            theme_noGrid()
        
    ggsave(paste0(out_prefix, tissue_1, "_", tissue_2,  "_mutation_strand_fc_scatter.pdf"), p, height = 5)
    
}

merge_tables <- function(x, y, min_muts = 6){
    
    merged <- rbind(x,y)
    merged <- merged %>%
        mutate(fc = log2(fc)) %>%
        filter(total >= min_muts) %>%
        select(gene, fc, mutation, tissue) %>%
        spread(tissue, fc)
    
    return(merged)
}

#perform_GO <- function(merged_fc, mut_type = "C>T", tissue_1, tissue_2, tissue1_greater_than = 1, tissue2_greater_than = 1) {
#    
#    merged_fc <- merged_fc[merged_fc$mutation == mut_type,]
#    
#    # getting the top/bottom fc genes
#    genes_test <- merged_fc[ ,tissue_1] >= tissue1_greater_than & merged_fc[ ,tissue_2] >= tissue2_greater_than
#    genes_test <- setNames(ifelse(genes_test, 0.001, 1), merged_fc$gene)
#    
#    # getting the top/bottom fc genes
#    genes_test_bottom <- merged_fc[ ,tissue_1] <= -tissue1_greater_than & merged_fc[ ,tissue_2] <= -tissue2_greater_than
#    genes_test_bottom <- setNames(ifelse(genes_test_bottom, 0.001, 1), merged_fc$gene)
#    
#    # GO enrichment
#    if(statistic == "ks") {
#        top_GO <- performGO(genes_test, statistic = "ks")
#        bottom_GO <- performGO(genes_test, gene_selection_function = function(p) p, statistic = statistic)
#    } else {
#        top_GO <- performGO(genes_test, gene_selection_function = function(p) p >= -log2(fdr_fc), statistic = statistic)
#        bottom_GO <- performGO(genes_test, gene_selection_function = function(p) p <= log2(fdr_fc), statistic = statistic)
#    }   
#    
#    # Get GO results
#    GO_tables_top <- get_GO_results(top_GO, test = statistic, fdr = F, bonf_signif = F, n = n_GO)
#    GO_tables_bottom <- get_GO_results(bottom_GO, test = statistic, fdr = F, bonf_signif = F, n = n_GO)
#    
#    # Get genes associated with GO
#    GO_genes_top <- get_genes_GO(genes = names(genes_test), x = GO_tables_top, GO_results = top_GO, n = n_GO)
#    GO_genes_bottom <- get_genes_GO(genes = names(genes_test), x = GO_tables_bottom, GO_results = bottom_GO, n = n_GO)
#    
#    GO_genes_top$gene_name <- queryRefCDS(x = GO_genes_top$gene, "gene_id", "gene_name")
#    GO_genes_bottom$gene_name <- queryRefCDS(x = GO_genes_bottom$gene, "gene_id", "gene_name")
#
#    
#    # Saving results
#    write.table(GO_tables_top, paste0(this_prefix, "_high_fc_GO.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
#    write.table(GO_tables_bottom, paste0(this_prefix, "_bottom_fc_GO.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
#    
#    write.table(GO_genes_top, paste0(this_prefix, "_high_fc_GO_genes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
#    write.table(GO_genes_top, paste0(this_prefix, "_bottom_fc_GO_genes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
#}



main()

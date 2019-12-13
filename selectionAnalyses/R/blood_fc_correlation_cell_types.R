# Correlates strand fold change values and cell type contribution
source("../../R/gtex.R", chdir = T)
source("../../R/plots.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
library("dplyr")
library("tidyr")
library("gplots")


main <- function(cmdArgs = commandArgs(T)) {
    
    out_prefix <- cmdArgs[1]
    fc_table_file <- cmdArgs[2]
    cibersort_file <- cmdArgs[3]
    
    #out_prefix <- "~/out_blood"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/n6_0.0_0.7/Whole_Blood-fcSrands.txt"
    #cibersort_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output"
    
    #out_prefix <- "~/out_blood"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/n6_0.0_0.7/Whole_Blood-fcSrands.txt"
    #cibersort_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Whole_Blood_all_genes.txt"
    
    #out_prefix <- "~/out_lung"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/n6_0.0_0.7/Lung-fcSrands.txt"
    #cibersort_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Lung_all_genes.txt"
    
    #out_prefix <- "~/out_lung_nk_cells"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung/n6_0.0_0.7/Lung-fcSrands.txt"
    #cibersort_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Lung_all_genes.txt"
    
    #out_prefix <- "~/out_lung_nk_cells_0.2"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung_0.2/n6_0.0_0.7/Lung-fcSrands.txt"
    #cibersort_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Lung_all_genes.txt"
    
    #out_prefix <- "~/out_lung_nk_cells_anti_0.2"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung_anti_0.2/n6_0.0_0.7/Lung-fcSrands.txt"
    #cibersort_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Lung_all_genes.txt"
    
    #out_prefix <- "~/out_lung_nk_cells_anti_0.1"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung_anti_0.1/n6_0.0_0.7/Lung-fcSrands.txt"
    #cibersort_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Lung_all_genes.txt"
    
    mutation_columns <- c("C_A", "C_T", "C_G", "T_A", "T_C", "T_G")
    
    # read files
    fc_table <- read.table(fc_table_file, sep = "\t", stringsAsFactors = F, header = T)
    cell_types <- read.table(cibersort_file, sep = "\t", stringsAsFactors = F, header = T)
    rownames(cell_types) <- gtexLongToSra(cell_types[,1])
    cell_types <- cell_types[cell_types$P.value < 0.05,]
    cell_types <- cell_types[,2:23]
    cell_type_columns <- colnames(cell_types)
    
    # append cell type to strand bias table
    fc_table <- cbind(fc_table, cell_types[fc_table$sample,])
    fc_table <- fc_table[,c(mutation_columns, cell_type_columns)]
    
    # Gathering data (tidy)
    fc_table <- fc_table %>%
        gather("mutation_type", "fc_strand", !!mutation_columns) %>%
        gather("cell_type", "percentage", !!cell_type_columns) %>%
        mutate(fc_strand = log2(fc_strand), percentage = 100* percentage)
    
    # Doing correlation
    cors <- fc_table %>%
        group_by(mutation_type, cell_type) %>%
        summarise(spearman = cor.test(fc_strand, percentage, method = "sp", na.rm = T)[["estimate"]],
                  pvalue = cor.test(fc_strand, percentage, method = "sp", na.rm = T)[["p.value"]]
                  ) %>%
        ungroup() %>%
        group_by(mutation_type) %>%
        mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
        filter(!is.na(spearman))
    
    cors <- cors[order(cors$mutation_type, cors$fdr),]
    
    
    # Getting matrix of correlation
    toPlot_heat <- get_matrix_heatmap(cors, column = spearman)
    pdf(paste0(out_prefix, "-heatmap_strand_bias_cell_types.pdf"))
    do_heatmap(toPlot_heat)
    dev.off()
    
    # Plotting the n highly correlated features as scatter plots
    plot_scatter(fc_table, cors, n_plots = 3, column = "spearman", x = "percentage", y = "fc_strand", out_prefix = out_prefix)
    
    # Saving correlation table
    write.table(cors, paste0(out_prefix, "-correlation_table_strand_bias_cell_types.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
    
        
}

get_matrix_heatmap <- function(cors, column) {
    
    column <- enquo(column)
    
    toPlot <- cors %>%
        dplyr::select(mutation_type, !!column, cell_type) %>%
        spread(mutation_type, !!column) %>%
        as.data.frame()
    
    rownames(toPlot) <- toPlot$cell_type
    toPlot <- as.matrix(toPlot[,-1])
    
    return(toPlot)
}

do_heatmap <- function(mat) {
    
    MUTATION_COLORS <- c(`C>A`="#3BC9F3", `C>G` = "#2B2E34", `C>T` = "#FC3218", `T>A` = "#CAD0CE", `T>C` = "#9CD169", `T>G` = "#F1CAC9")
    
    .scale <- "none"
    .col <- colorRampPalette(c("blue", "black", "yellow"))(13)
    
    # Plot heatmap
    gplots::heatmap.2(mat, Rowv = T, Colv = F, trace = "none", ColSideColors = MUTATION_COLORS, col = .col, scale = .scale, density.info = "none")
    
}

plot_scatter <- function(fc_table, cors, n_plots = 3, column = "spearman", x = "percentage", y = "fc_strand", out_prefix) {
    for(mut_type in unique(cors$mutation_type)) {
        
        current_cors <- cors[cors$mutation_type == mut_type,]
        current_cors <- current_cors[order(current_cors[,column, drop = T]),]
        
        #Select top features
        cell_types <- c("NK.cells.resting", "T.cells.CD8", "Neutrophils", head(current_cors$cell_type, n_plots)[n_plots], tail(current_cors$cell_type, n_plots))
        cell_types <- cell_types[1:(n_plots*2)]
        
        current_fc <- fc_table[fc_table$mutation_type == mut_type & fc_table$cell_type %in% cell_types,]
        
        p <- scatter(current_fc, x = x, y = y, facet_y = "cell_type", ncolFactor = 2, regression = T, method = "sp", 
                     alpha = 0.5, xlab = "Cell type content", ylab = paste0(mut_type, " log2(coding/non-coding)"),
                     pSize = 0.5) +
              theme_noGrid()
        
        ggsave(paste0(out_prefix, "-scatter_strand_bias_cell_types_", mut_type, ".pdf"), p, height = 6, width = 4)
    }

}

main()

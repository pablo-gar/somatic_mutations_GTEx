# Correlates strand fold change values and neutrophil expression in lung
# 1. Get spearman correlation for each sample between their expression profile and neutrophil expression
# 2. Correlata the above correlation with strand bias
source("../../R/gtex.R", chdir = T)
source("../../R/plots.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
source("../../R/geneTools.R", chdir = T)
library("dplyr")
library("tidyr")
library("gplots")


main <- function(cmdArgs = commandArgs(T)) {
    
    out_prefix <- cmdArgs[1]
    fc_table_file <- cmdArgs[2]
    cibersort_LMM22_file <- cmdArgs[3]
    
    #out_prefix <- "out_"
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/n6_0.0_0.7/Lung-fcSrands.txt"
    #cibersort_LMM22_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/LM22.txt.ready"
    
    mutation_columns <- c("C_A", "C_T", "C_G", "T_A", "T_C", "T_G")
    
    # read files
    fc_table <- read.table(fc_table_file, sep = "\t", stringsAsFactors = F, header = T)
    cell_types <- read.table(cibersort_LMM22_file, sep = "\t", stringsAsFactors = F, header = T)
    
    # Read gtex expression
    fc_table$gtexInd <- sraToGtex(fc_table$sample, formatOut = "long")
    exp_mat <- readExpression(fc_table$gtexInd)
    
    # Select gene expression of neutrophil genes
    cell_types <- cell_types[cell_types$Gene.symbol %in% rownames(exp_mat),]
    exp_mat <- exp_mat[cell_types$Gene.symbol,]
    
    # Get correlations per individual
    all_cors <- setNames(rep(0, ncol(cell_types) - 1), colnames(cell_types)[-1])
    for(i in names(all_cors)) {
        cell_type_cor <- apply(exp_mat, 2, function(x) cor(x, cell_types[,i], method = "sp"))
        all_cors[i] <- cor(cell_type_cor, fc_table$C_T, method = "sp")
        
        toPlot <- data.frame(x = cell_type_cor, y = fc_table$C_T)
        p <- scatter(toPlot, x = "x", y = "y",  method_cor = "sp", regression = T, pSize = 1) + theme_noGrid()
        ggsave(paste0(out_prefix, "scatter_lung_C_T_strand_asymmetry_", i, ".pdf"), p, width = 4, height = 4)
    }
    
    all_cors <- sort(all_cors)
    write.table(as.data.frame(all_cors), paste0(out_prefix, "C_T_strand_asymmetry_cors.txt"), sep = "\t", col.names = F, row.names = T, quote = F)
}

readExpression <- function(gtexInd) {
    
    exp_mat <- readAllGtexExpression(gtexInd)
    exp_mat[,1] <- gsub("\\.\\d+","",exp_mat[,1])
    hugo_symbols <- ensemblToAlias(exp_mat[,1])
    exp_mat[,1] <- hugo_symbols
    exp_mat <- exp_mat[!is.na(exp_mat)[,1],]
    
    exp_mat <- exp_mat[!duplicated(exp_mat[,1]),]
    
    rownames(exp_mat) <- exp_mat[,1]
    exp_mat <- exp_mat[,-1]
    
    return(exp_mat)
    
}

main()

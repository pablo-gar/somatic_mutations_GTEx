source("../../R/plots.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
source("../../R/misc.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    mut_file <- cmdArgs[1]
    fc_strand_file <- cmdArgs[2]
    out_prefix <- cmdArgs[3]
    
    mutation_types <-  c("C_A", "C_G", "C_T", "T_A", "T_C", "T_G")
    
    artifacts <- c("SMRIN", "SMTSISCH", "SMTSPAX",  "n_uniqueMapped", "transcriptome_diversity")
    
    #mut_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Lung-n6_0.0_0.7.txt"
    #fc_strand_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/n6_0.0_0.7/Lung-fcSrands_allInds_noMetadata.txt"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/strand_fc/mutationSignature_associations_mutationLoad/Lung/n6_0.0_0.7/"
    
    mut_covariates <- read.table(mut_file, sep = "\t", stringsAsFactors = F, header = T)
    mut_covariates <- get_residuals(mut_covariates, artifacts, mutation_types)
    fc_strand <- read.table(fc_strand_file, sep = "\t", stringsAsFactors = F, header = T)
    
    
    # Perform correlations
    results <- list()
    counter <- 1
    for(mut_type in mutation_types) {
        
        current_mutation_load <- mut_covariates[mut_covariates$sample %in% fc_strand$sample, ]
        current_fc_strand <- setNames(fc_strand[ fc_strand$sample %in% current_mutation_load$sample, mut_type], fc_strand$sample[fc_strand$sample %in% current_mutation_load$sample])
        current_fc_strand <- log2(current_fc_strand)
        
        cor_results <- cor.test(current_mutation_load[,mut_type], current_fc_strand[current_mutation_load$sample], method = "sp")
        
        toPlot <- data.frame(mutation_load = current_mutation_load[,mut_type], strand_fc = current_fc_strand)
        p <- scatter(toPlot, y = "strand_fc", x = "mutation_load", method_cor = "sp", regression = T)
        p <- p + 
            xlab(paste0("Normalized mutations (", mut_type, ")")) +
            ylab(paste0(mut_type, " strand bias\nlog2(transcribed/non−transcribed)")) +
            theme_noGrid()
        
        ggsave(paste0(out_prefix, "mut_type_", mut_type, "_scatter.pdf"), p, height = 4, width = 4)
        
        results[[counter]] <- data.frame(mut_type = mut_type, spearman = cor_results$est, pvalue = cor_results$p.value, stringsAsFactors = F)
        counter <- counter + 1
    
    }
    
    results <- do.call(rbind, results)
    results <- results[order(results$mut_type, results$pvalue),]
    
    write.table(results, paste0(out_prefix, "1_stats.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
    
}

get_residuals <- function(x, artifacts, mut_columns) {
    
    for(mut in mut_columns) {
        x$temp <- 0
        x[,mut] <- rankitNormalize_vector(x[,mut])
        x <- x[rowSums(is.na(x)) == 0,]
        x[,mut] <-  resid(lm(as.formula(paste0(mut, "~.")), data = x[,colnames(x) %in% c(mut, artifacts)]))
    }
    
    return(x)
}


main()

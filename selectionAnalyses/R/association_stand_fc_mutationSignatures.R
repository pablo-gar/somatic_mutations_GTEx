source("../../R/plots.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    fc_table_file <- cmdArgs[1]
    signature_weight_file <- cmdArgs[2]
    out_prefix <- cmdArgs[3]
    percentage <- as.logical(cmdArgs[4]) # Calculate percentage across mutations ?
    
    mutation_types <-  c("C_A", "C_G", "C_T", "T_A", "T_C", "T_G")
    
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/n6_0.0_0.7/Lung-fcSrands_allInds_noMetadata.txt"
    #signature_weight_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/reverse_cancer_signatures_cors/Lung-n6_0.0_0.7_H_cancer_signatures.txt"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/strand_fc/mutationSignature_associations_cor/Lung/n6_0.0_0.7/"
    #percentage <- F
    
    #fc_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/n6_0.0_0.7/Lung-fcSrands_allInds_noMetadata_total_mutations.txt"
    #signature_weight_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/reverse_cancer_signatures_cors/Lung-n6_0.0_0.7_H_cancer_signatures.txt"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/strand_fc/mutationSignature_associations_cor_totalMuts/Lung/n6_0.0_0.7/"
    #percentage <- T
    
    fc_table <- read.table(fc_table_file, sep = "\t", stringsAsFactors = F, header = T)
    signature_weight <- read.table(signature_weight_file, sep = "\t", stringsAsFactors = F, header = T)
    
    if(percentage) {
        fc_table[,mutation_types] <- fc_table[,mutation_types] / colSums(fc_table[,mutation_types], na.rm = T)
    }
    
    
    # Perform correlations
    results <- list()
    counter <- 1
    for(mut_type in mutation_types) {
        for (mut_sign in signature_weight$signature) {
            
            current_fc <- fc_table[fc_table$sample %in% colnames(signature_weight), ]
            
            if(!percentage)
                current_fc[,mut_type] <- log2(current_fc[,mut_type])
            
            current_sign <- as.numeric(signature_weight[ signature_weight$signature == mut_sign, current_fc$sample])
            cor_results <- cor.test(current_fc[,mut_type], current_sign, method = "sp")
            
            toPlot <- data.frame(strand_fc = current_fc[,mut_type], signature_w = current_sign)
            p <- scatter(toPlot, y = "strand_fc", x = "signature_w", method_cor = "sp", regression = T)
            p <- p + 
                xlab(paste0("Signature ", mut_sign, " weight")) +
                ylab(paste0(mut_type, " strand bias\nlog2(transcribed/non−transcribed)")) +
                theme_noGrid()
            
            ggsave(paste0(out_prefix, "mut_type_", mut_type, "_signature_", mut_sign, "_scatter.pdf"), p, height = 4, width = 4)
            
            results[[counter]] <- data.frame(mut_type = mut_type, signature = mut_sign, spearman = cor_results$est, pvalue = cor_results$p.value, stringsAsFactors = F)
            counter <- counter + 1
        }
    }
    
    results <- do.call(rbind, results)
    results <- results[order(results$mut_type, results$pvalue),]
    
    write.table(results, paste0(out_prefix, "1_stats.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
    
}

main()

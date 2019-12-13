j Reads multiple files from get_mutation_residual_afterArtifactRemoval.R and 
# compares mutations according to certain covariates
# Rscript compare_two_tissues_mutations.R /scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/race_skin_comparisons/out_ /scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Skin_Sun_Exposed_Lower_leg-n6_0.0_0.7.txt /scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Skin_Not_Sun_Exposed_Suprapubic-n6_0.0_0.7.txt

library("ggplot2")
library("dplyr")
library("tidyr")
source("../../R/gtex.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)

covariates_compare <- "RACE"
artifacts <- c("SMRIN", "SMTSISCH", "SMTSPAX",  "n_uniqueMapped", "transcriptome_diversity")
race <- setNames(c("Asian", "African American", "Caucasian", "Native American", "Not Reported", "Unknown"),c(1,2,3,4,98,99))

main <- function(cmdArgs = commandArgs(T)) {
    
    outPrefix <- cmdArgs[1]
    mut_files <- cmdArgs[2:length(cmdArgs)]
    
    #mut_files <- file.path("/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates", 
    #                       c("Skin_Sun_Exposed_Lower_leg-n6_0.0_0.7.txt", "Skin_Not_Sun_Exposed_Suprapubic-n6_0.0_0.7.txt")
    #                       )
    #outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/skin_caucasian_african/n6_0.0_0.7_"
    
    mut_covariates <- read_mutFiles(mut_files)
    mut_covariates <- mut_covariates[rowSums(is.na(mut_covariates)) ==0,]
    mut_covariates$gtex_id <- sraToGtex(mut_covariates$sample, formatOut = "short")
    mut_covariates$RACE <- readMetadataGTEX("RACE")[mut_covariates$gtex_id,]
    mut_covariates <- mut_covariates[mut_covariates$RACE %in% c(2,3),]
    mut_covariates$RACE <- race[as.character(mut_covariates$RACE)]
    
    # Gets residuals on C>T mutations
    mut_covariates$C_T_resid <- 0
    for(tissue in unique(mut_covariates$tissue)) {
        current <- mut_covariates[mut_covariates$tissue == tissue,]
        current$C_T_resid <- resid(lm(C_T~., data = current[,c("C_T", artifacts)]))
        mut_covariates[mut_covariates$tissue == tissue,] <- current
    }
    
    # Plots C>T residuals in different tissues and races
    p_resid <- ggplot(mut_covariates, aes(x = tissue, y = C_T_resid)) +
        geom_boxplot(aes(fill = RACE)) +
        theme_grid_y()
        
    
    # Calculates freqs
    mut_covariates[,1:9] <- mut_covariates[,1:9] / mut_covariates[,9] * 100
    test_results <- mut_covariates %>%
        group_by(RACE) %>%
        summarise(label = paste0("p=", signif(wilcox.test(C_T~tissue)[["p.value"]], 2)), C_T = max(C_T) + 1) %>%
        ungroup()
    
    p_freq_tissue_race <- ggplot(mut_covariates, aes(x = RACE, y = C_T)) +
                             geom_violin(aes(fill = tissue), width = 0.30) +
                             geom_boxplot(aes(fill = tissue), width = 0.25) +
                             scale_fill_manual(values = c("#d8be88", "#ffbb30")) +
                             geom_text(aes(label = label), data = test_results, fontface = "italic") +
                             xlab("") + ylab("C>T mutations (%)") +
                             theme_grid_y() +
                             theme(legend.position = "top")
    
    # Calculate fold-change, select only the ones for which we have both tissues
    fc_mutations <- mut_covariates %>%
        select(C_T, RACE, gtex_id, tissue) %>%
        spread(tissue, C_T)
    fc_mutations <- fc_mutations[rowSums(is.na(fc_mutations)) == 0,]
    fc_mutations$fc <- log2(fc_mutations[,ncol(fc_mutations)] /  fc_mutations[,ncol(fc_mutations)-1])
    test_results <- paste0("p=", signif(wilcox.test(fc~RACE, data = fc_mutations)[["p.value"]], 2))
    
    p_foldChange <- ggplot(fc_mutations, aes(x = RACE, y = fc)) +
                       geom_violin(fill = "grey80", width = 0.5) +
                       geom_boxplot(fill = "grey75", width = 0.25) +
                       annotate("text", label = test_results, x = -Inf, y = Inf, hjust = 0, vjust = 1, fontface ="italic") +
                       xlab("") + ylab("C>T Fold Change [log2(Sun/Non)]") +
                       theme_grid_y() +
                       theme(axis.text.x = element_text(angle = 20, hjust = 1))
                   
    #saving results
    ggsave(paste0(outPrefix, "C_T_foldChange_by_race.pdf"), p_foldChange, width = 2, height = 5.5)
    ggsave(paste0(outPrefix, "C_T_frequency_by_tissue_and_race.pdf"), p_freq_tissue_race, width = 2.7, height = 5.5)

    
}

read_mutFiles <- function(mut_files) {
    results <- list()
    for(i in mut_files) {
        results[[i]] <- read.table(i, sep = "\t", stringsAsFactors = F, header = T)
        results[[i]]$tissue <- gsub(".txt", "", basename(i))
        results[[i]]$tissue <- gsub("-.+", "", results[[i]]$tissue)
    }
    results <- do.call(rbind, results)
    rownames(results) <- 1:nrow(results)
    return(results)
}

main()

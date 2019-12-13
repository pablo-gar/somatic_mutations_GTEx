# Reads multiple files from get_mutation_residual_afterArtifactRemoval.R and 
# compares mutations according to certain covariates
# Rscript compare_two_tissues_mutations.R /scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/race_skin_comparisons/out_ /scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Skin_Sun_Exposed_Lower_leg-n6_0.0_0.7.txt /scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Skin_Not_Sun_Exposed_Suprapubic-n6_0.0_0.7.txt

library("ggplot2")
library("dplyr")
library("tidyr")
source("../../R/gtex.R", chdir = T)
source("../../R/plots.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
source("../../R/misc.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    

    outPrefix <- cmdArgs[1]
    variables <- cmdArgs[2]
    mut_files <- cmdArgs[3:length(cmdArgs)]
    
    artifacts <- c("SMRIN", "SMTSISCH", "SMTSPAX",  "n_uniqueMapped", "transcriptome_diversity")
    race <- setNames(c("Asian", "African American", "Caucasian", "Native American", "Not Reported", "Unknown"),c(1,2,3,4,98,99))
    MUT_COLUMNS <- 1:7
    
    #mut_files <- file.path("/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates", 
    #                       c("Spleen-n6_0.0_0.7.txt", "Pancreas-n6_0.0_0.7.txt")
    #                       )
    #variables <- c("BMI")
    #outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutationsVsMetadata_scatter_box/plots_"
    
    #mut_files <- file.path("/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates", 
    #                       c("Breast_Mammary_Tissue-n6_0.0_0.7.txt")
    #                       )
    #variables <- c("GENDER")
    #outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutationsVsMetadata_scatter_box/n6_0.0_0.7/Breast_gender/plots_"
    
    #mut_files <- file.path("/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates", 
    #                       c("Spleen-n6_0.0_0.7.txt", "Whole_Blood-n6_0.0_0.7.txt", "Skin_Sun_Exposed_Lower_leg-n6_0.0_0.7.txt", "Adrenal_Gland-n6_0.0_0.7.txt", "Brain_Hypothalamus-n6_0.0_0.7.txt", "Small_Intestine_Terminal_Ileum-n6_0.0_0.7.txt")
    #                       )
    
    #mut_files <- list.files("/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates", full = T)
    #mut_files <- mut_files[!grepl("EXO", mut_files)]
    #variables <- c("AGE", "GENDER", "BMI")
    #outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutationsVsMetadata_scatter_box/plots_"
    
    mut_covariates <- read_mutFiles(mut_files)
    
    # Gets residuals on all mutations
    mut_covariates <- get_residuals(mut_covariates, artifacts, MUT_COLUMNS)
    
    # Plots
    save_plots(mut_covariates, MUT_COLUMNS, variables, outPrefix)
    
}

read_mutFiles <- function(mut_files, ignore_gender = F) {
    results <- list()
    for(i in mut_files) {
        flush.console()
        print(i)
        results[[i]] <- read.table(i, sep = "\t", stringsAsFactors = F, header = T)
        results[[i]]$tissue <- gsub(".txt", "", basename(i))
        results[[i]]$tissue <- gsub("-.+", "", results[[i]]$tissue)
        results[[i]] <- results[[i]] [,!colnames(results[[i]]) %in% c("C_C", "T_T")]
    }
    
    #shared_colnames <- Reduce(intersect, lapply(results, colnames))
    #results <- lapply(results, function(x, col.names) x[,col.names], col.names = shared_colnames)
    
    #results <- do.call(rbind, results)
    #rownames(results) <- 1:nrow(results)
    return(results)
}

# Get residuals per tissue per mutations type
get_residuals <- function(y, artifacts, mut_columns) {
    
    for(j in seq_along(y)) {
        x <- y[[j]]
        for(i in mut_columns) {
            mut <- colnames(x)[i]
            x$temp <- 0
            x[,mut] <- rankitNormalize_vector(x[,mut])
            x <- x[rowSums(is.na(x)) == 0,]
            x$temp <-  resid(lm(as.formula(paste0(mut, "~.")), data = x[,colnames(x) %in% c(mut, artifacts)]))
            colnames(x)[ncol(x)] <- paste0(mut, ".resid")
        }
        y[[j]] <- x
    }
    
    return(y)
}

save_plots <- function(x, mut_columns, variables, outPrefix) {
    
    
    results <- list()
    counter <- 1
    
    for(j in seq_along(x)) {
        toPlot <- x[[j]]
        tissue <- toPlot$tissue[1]
        mut_resid_columns <- (ncol(toPlot) - length(mut_columns) + 1) : ncol(toPlot)
        
        for(i in mut_resid_columns) {
            mut <- colnames(toPlot)[i]
        
            for(variable in variables) {
                # Separating in two groups for boxplots
                if (variable == "GENDER") {
                    toPlot$group <- as.factor(toPlot[,variable])
                    toPlot_final <- toPlot
                } else {
                    nBins <- 4
                    bins <- seq(0,1, length.out = nBins + 1)
                    
                    quants <- quantile(toPlot[,variable], bins)
                    label <-  paste0(quants[-(nBins+1)], "-", quants[-1])
                    groups <- as.character(cut(toPlot[,variable],  quants, label, include.lowest = T))
                    toPlot$group <- groups
                    toPlot_final <- toPlot[toPlot$group %in% label[c(1,nBins)],]
                }
                
                if (length(unique(toPlot_final$group)) != 2) {
                    cat(paste("There are not 2 groups for", variable, tissue, mut, "\n"))
                    next
                }
                
                # Plotting
                spearman_results <- cor.test(toPlot[,variable], toPlot[,mut], method = "sp")
                p_scatter <- scatter(toPlot, x = variable, y = mut, method_cor = "sp", regression = T) + 
                    theme_noGrid() +
                     xlab(variable) +
                    ylab(paste0("Normalized ", mut, "mutations")) 
                
                wilcox_test_resuls <- signif(wilcox.test(as.formula(paste0(mut, "~group")), data = toPlot_final)[["p.value"]], 2)
                colour_fill <- colorRampPalette(c("#657691", "#b5d1ff"))(length(unique(toPlot_final$group)))
                p_box <- ggplot(toPlot_final, aes_string(y = mut, x = "group")) +
                     geom_violin(fill="grey85", width = 0.6)+
                     geom_boxplot(fill = colour_fill, width = 0.25)+
                     ylab(paste0("Normalized ", mut, "mutations")) +
                     xlab(variable) +
                     annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
                                  label = paste0("p=", signif(wilcox_test_resuls)), fontface = "italic") +
                    theme_grid_y() 
                
                long_prefix <- paste0(outPrefix, tissue, "_", variable, "_", mut, "_")
                ggsave(paste0(long_prefix, "scatter.pdf"), p_scatter, height = 4, width = 4)
                ggsave(paste0(long_prefix, "box.pdf"), p_box, height = 5, width = 2)
                
                
                # Saving numeric results
                results[[counter]] <- data.frame(tissue = tissue, mut_type = mut, feature = variable,
                                                 spearman.rho = spearman_results[["estimate"]], spearman.p = spearman_results[["p.value"]],
                                                 wilcox.p = wilcox_test_resuls)
                counter <- counter + 1
                
            }
        }
    }
    
    results <- do.call(rbind, results)
    results <- results %>%
        group_by(feature, mut_type) %>%
        mutate(wilcox.fdr_bh = p.adjust(wilcox.p, method = "BH"), spearman.fdr_bh = p.adjust(wilcox.p, method = "BH")) %>%
        ungroup()
    
    write.table(results, paste0(outPrefix, "stats.txt"), sep = "\t", quote = F, row.names = F)
    
}
main()

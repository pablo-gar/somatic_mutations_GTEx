library("ggplot2")
library("dplyr")
library("tidyr")
source("../../R/ggthemes.R")
source("../../R/plots.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    workingDir <- cmdArgs[1] # Dir with counts, e.g mutationCount/count/
    maf <- cmdArgs[2] # Maf of interest 
    out_prefix <- cmdArgs[3]
    sample_table <- cmdArgs[4]
    
    metadata_cols <- c("AGE", "BMI", "GENDER", "SMRIN", "SMTSISCH", "SMTSPAX", "C1", "C2", "C3", "n_uniqueMapped")
    ignore_metadata_cols <- c("SMRIN", "SMTSISCH", "SMTSPAX",  "C1", "C2", "C3", "RACE", "transcriptome_diversity") 
    metadata_cols_interest <- metadata_cols[!metadata_cols %in% ignore_metadata_cols]
    
    #workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/"
    #maf <- "n6_0.0_0.7"
    #out_prefix <- "/home/users/paedugar/deleteme_"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/regressions_seq_depth/"
    #sample_table <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/sample_table.txt"
    
    
    # Read giles
    samples <- read.table(sample_table, sep = "\t", header = T, stringsAsFactors = F)
    samples <- samples[,!colnames(samples) %in% ignore_metadata_cols]
    
    # Get mutatqon count per sample 
    mutation_counts <- read_counts_samples(samples, workingDir, maf)
    rownames(mutation_counts) <- mutation_counts$sraId
        
    # Getting linear models for each mutations type: mutations ~ read_depth
    mutation_models <- get_models(mutation_counts, metadata_cols_interest)
    
    # Doing plots per tissue
    all_muts <- list()
    for(i in names(mutation_models)) {
        # Getting residuals and predicted values
        resids <- resid(mutation_models[[i]])
        predicted <- predict(mutation_models[[i]])
        
        # Adding to data frame and reshaping
        current <- mutation_counts[names(resids), ]
        current$residual <- resids
        current$predicted <- predicted
        
        current <- gather(current, "variable", "value", residual, predicted, n_uniqueMapped)
        
        p <- scatter(current, x = "mutations", y = "value", facet_x = "variable", labelSize = 0, pSize = 2, alpha = 0.5) + 
            theme_noGrid()
        ggsave(paste0(out_prefix, i, "_scatter_seqDepth_predicted_resid.pdf"), p, width = 8, height = 3)
        
        all_muts[[i]] <- current
    }
            
    all_muts <- do.call(rbind, all_muts)
    
    
    toEliminate <- all_muts[ all_muts$variable == "residual" & all_muts$value >= 1500, "sraId"]
    all_muts$pass <- apply(all_muts, 1, function(x) !x["sraId"] %in% toEliminate)
        
    p <- scatter(all_muts, x = "mutations", y = "value", pColour = "pass", facet_x = "variable", labelSize = 0, pSize = 2, alpha = 0.1) + 
        theme_noGrid()
    
    
    all_muts <- spread(all_muts, variable, value)

    # Save results
    ggsave(paste0(out_prefix, "1_allTissues_scatter_seqDepth_predicted_resid.pdf"), p, width = 8, height = 3)
    write.table(all_muts, paste0(out_prefix, "1_allTissues_stats.txt"), sep = "\t", quote = F)
    
    # Print samples to be removed
    cat(paste0(all_muts$sraId[!all_muts$pass], collapse = " "))

        
}

read_counts_samples <- function(samples, workingDir, maf) {
    
    results <- list()
    counter <- 1
    for(i in 1:nrow(samples)) {
        current_file <- file.path(workingDir, samples$tissue[i], maf, paste0(samples$sraId[i], ".txt"))
        
        if(!file.exists(current_file))
           next 
           
        counts <- read.table(current_file, header = F, sep = "\t")
        current <- data.frame(
                              mutations = sum(counts), 
                              C_A = counts[2,1], C_T = counts[2,2], C_G = counts[2,3],
                              T_A = counts[1,1], T_C = counts[1,4], T_G = counts[1,3]
                              )
        
        rownames(current) <- samples$sraId[i]
        
        results[[counter]] <- current
        counter <- counter + 1 
    }
    
    results <- do.call(rbind, results)
    
    samples <- samples[samples$sraId %in% rownames(results),]
    samples <- cbind(samples, results[samples$sraId,])
    
    return(samples)
}

get_models <- function(mutation_counts, columns = c("AGE", "BMI", "GENDER", "SMRIN", "SMTSISCH", "SMTSPAX", "C1", "C2", "C3", "n_uniqueMapped", "transcriptome_diversity")) {
    
    models <- list()
    for(i in unique(mutation_counts$tissue)) {
        ids <- mutation_counts[mutation_counts$tissue == i, "sraId"]
        current <- mutation_counts[mutation_counts$tissue == i, c("mutations", columns)]
        rownames(current) <- ids
        form <- as.formula(paste("mutations ~ ."))
        models[[i]] <- lm(form, data = current, na.action = na.exclude)
    }
    
    return(models)
    
}

main()

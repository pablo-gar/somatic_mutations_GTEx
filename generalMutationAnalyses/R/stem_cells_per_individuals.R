library("ggplot2")
library("dplyr")
library("ggrepel")
library("metap")
source("../../R/ggthemes.R")
source("../../R/misc.R")
source("../../R/FDR.R")
source("../../R/plots.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    workingDir <- cmdArgs[1] # Dir with counts, e.g mutationCount/count/
    maf <- cmdArgs[2] # Maf of interest 
    out_prefix <- cmdArgs[3]
    n_tissues_cutoff <- as.numeric(cmdArgs[4])
    stem_cell_division_file <- cmdArgs[5]
    sample_table <- cmdArgs[6]
    
    metadata_cols <- c("AGE", "BMI", "GENDER", "SMRIN", "SMTSISCH", "SMTSPAX", "C1", "C2", "C3", "n_uniqueMapped", "transcriptome_diversity")
    ignore_metadata_cols <- c("SMRIN", "SMTSISCH", "SMTSPAX",  "C1", "C2", "C3", "RACE") 
    metadata_cols_interest <- metadata_cols[!metadata_cols %in% ignore_metadata_cols]
    tissues_to_remove <- c("Esophagus_Muscularis")
    
    #workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/"
    #maf <- "n6_0.0_0.7"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/totalCountsBox/raw_n6_0.0_0.7_"
    #n_tissues_cutoff <- 4
    #stem_cell_division_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/stem_cell_divisions.txt"
    #sample_table <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/sample_table.txt"
    
    
    # Read giles
    stem_cell <- read_steam_cell_file(stem_cell_division_file)
    samples <- read.table(sample_table, sep = "\t", header = T, stringsAsFactors = F)
    samples <- samples[samples$RACE == 3, ] # only caucasians
    samples <- samples[,!colnames(samples) %in% ignore_metadata_cols]
    samples <- samples[!samples$tissue %in% tissues_to_remove,]
    
    
    # Get mutation count per sample 
    mutation_counts <- read_counts_samples(samples, workingDir, maf)
    
    # Get mutation count per  tissue
    mutation_counts_tissue <- get_counts_per_tissue(mutation_counts)
    
    # Removing samples with NAs
    mutation_counts <- mutation_counts[rowSums(is.na(mutation_counts)) == 0,]
        
    # Getting linear models for each mutations type: mutations ~ read_depth
    mutation_models_tissue <- get_models_tissue(mutation_counts_tissue, "n_uniqueMapped")
    mutation_models <- get_models(mutation_counts, metadata_cols_interest)
    
        
    # Getting individuals with at least n tissues
    tissue_count <- get_individuals_most_tissues_stem(mutation_counts, stem_cell, merge_patterns = c("Brain"), to = c("Brain_Cortex"))
    gtexInds <- tissue_count[tissue_count$tissues_stem >= n_tissues_cutoff, "gtexInd", drop = T]
    
    mutation_counts <- mutation_counts[mutation_counts$gtexInd %in% gtexInds, ]
    mutation_counts <- mutation_counts[mutation_counts$tissue %in% rownames(stem_cell), ]
    
    # Append stem cell info
    mutation_counts$log10_stem_cell_division <- stem_cell[mutation_counts$tissue, "log10_stem_cell_division"]
    
    # Working in all mutation types
    all_cors <- list()
    fdrs <- list()
    for(i in c("mutations", "T_A", "T_G", "T_C", "C_A", "C_T", "C_G")) {
        
        # Getting obs/predicted mutations
        predicted_muts <- predict(mutation_models[[i]], mutation_counts) 
        mutation_counts$log2_fc_obs_expected <- log2( (mutation_counts$mutations + abs(min(predicted_muts)) + 1) / (predicted_muts + abs(min(predicted_muts)) + 1))
        
        # Getting averages for repeated tissues
        mutation_counts_averaged <- get_average(mutation_counts, column_ind = "gtexInd", tissue_pattern = "Brain", new_name = "Brain", column_average = "log2_fc_obs_expected", 
                                                other_cols = c("gtexInd", "log10_stem_cell_division", "AGE", "tissue"))
        
        spear_cors <- mutation_counts_averaged %>%
            group_by(gtexInd) %>%
            summarise(spearman = cor.test(log10_stem_cell_division, log2_fc_obs_expected, method = "sp")[["estimate"]],
                      pvalue = cor.test(log10_stem_cell_division, log2_fc_obs_expected, method = "sp")[["p.value"]],
                      n_obs = n(),
                      age = AGE[1]) %>%
            ungroup()
        
        fdr <- cor_by_group_fdr(mutation_counts_averaged, x = log10_stem_cell_division, y = log2_fc_obs_expected, group_col = gtexInd, 1000, 0.05, theoretical_fdr = F)
        
        spear_cors$mut_type <-i 
        all_cors[[i]] <- spear_cors
        fdrs[[i]] <- data.frame(fdr = fdr, mut_type = i, stringsAsFactors = F)
    }
    
    # Saving results
    all_cors <- do.call(rbind, all_cors)
    fdrs <- do.call(rbind, fdrs)
    all_cors <- all_cors[order(all_cors$mut_type, -all_cors$spearman),]
    write.table(all_cors, paste0(out_prefix, "_stem_cell_cors.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
    write.table(fdrs, paste0(out_prefix, "_stem_cell_cors_fdrs.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
    
    
    ###########################################
    # save plot for a good individuals
    ind <- "GTEX-ZUA1"
    mut_type <- "C_T"
    
    # Getting obs/predicted mutations
    predicted_muts <- predict(mutation_models[[mut_type]], mutation_counts) 
    mutation_counts$log2_fc_obs_expected <- log2( (mutation_counts$mutations + abs(min(predicted_muts)) + 1) / (predicted_muts + abs(min(predicted_muts)) + 1))
    
    # Getting averages for repeated tissues
    mutation_counts_averaged <- get_average(mutation_counts, column_ind = "gtexInd", tissue_pattern = "Brain", new_name = "Brain", column_average = "log2_fc_obs_expected", 
                                            other_cols = c("gtexInd", "log10_stem_cell_division", "AGE", "tissue"))
    
    toPlot <- mutation_counts_averaged[mutation_counts_averaged$gtexInd == ind,]
    
    p <- scatter(toPlot, x = "log10_stem_cell_division", y = "log2_fc_obs_expected", method = "sp", regression = T) + 
        geom_text_repel(aes(label = tissue)) +
        xlab("log10(stem cell divisions)") +
        ylab("Mutations log2(observed/expected)") +
        ggtitle(ind) +
        theme_noGrid()
    ggsave(paste0(out_prefix, ind, "_scatter.pdf"), p, width = 4, height = 4)

    ###############################
    

        
}

read_steam_cell_file <- function(stem_cell_division_file) {
    
    # Read stem cell division per sample
    stem_cell <- read.table(stem_cell_division_file, sep = "\t", header = T, stringsAsFactors = F)
    stem_cell <- stem_cell[!is.na(stem_cell$cancer_type),]
    stem_cell$log10_stem_cell_division <- log10(stem_cell$n_divisions_all_stem_cells_per_lifetime)
    stem_cell$log10_life_time_cancer_risk <- log10(stem_cell$life_time_cancer_risk)
    
    rownames(stem_cell) <- stem_cell$tissue
    
    return(stem_cell)
}

get_individuals_most_tissues_stem <- function(samples, stem_cell, merge_patterns = c("Brain"), to = c("Brain_Cortex")) {
    
    for(i in seq_along(merge_patterns)) {
        samples[grepl(merge_patterns[i], samples$tissue), "tissue"] <- to[i]
    }
    
    stem_tissue_count <- samples %>%
        group_by(gtexInd) %>%
        summarise(tissues_stem = sum(unique(tissue) %in% stem_cell$tissue)) %>%
        ungroup()
        
    stem_tissue_count <- stem_tissue_count[order(-stem_tissue_count$tissues_stem),]
    return(stem_tissue_count)
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

get_counts_per_tissue <- function(mutation_counts) {
    
    mutation_counts_tissue <- mutation_counts %>%
        group_by(tissue) %>%
        summarise(n_uniqueMapped = sum(as.numeric(n_uniqueMapped)), mutations = sum(as.numeric(mutations)), 
                  C_A = sum(as.numeric(C_A)), C_T = sum(as.numeric(C_T)), C_G = sum(as.numeric(C_G)), 
                  T_A = sum(as.numeric(T_A)), T_C = sum(as.numeric(T_C)), T_G = sum(as.numeric(T_G)),
                  ) %>%
        ungroup()
    
    return(mutation_counts_tissue)
        
    
}

get_models_tissue <- function(mutation_counts_tissue, read_count_column) {
    
    mutation_counts_tissue <- as.data.frame(mutation_counts_tissue)
    models <- list()
    for(i in c("mutations", "T_A", "T_G", "T_C", "C_A", "C_T", "C_G")) {
        form <- as.formula(paste(i, "~", read_count_column))
        models[[i]] <- lm(form, data = mutation_counts_tissue)
        
    }
    
    return(models)
    
}

get_models <- function(mutation_counts, columns = c("AGE", "BMI", "GENDER", "SMRIN", "SMTSISCH", "SMTSPAX", "C1", "C2", "C3", "n_uniqueMapped", "transcriptome_diversity")) {
    
    models <- list()
    for(i in c("mutations", "T_A", "T_G", "T_C", "C_A", "C_T", "C_G")) {
        current <- mutation_counts[,c(i, columns)]
        form <- as.formula(paste(i, "~ ."))
        models[[i]] <- lm(form, data = current, na.action = na.exclude)
    }
    
    return(models)
    
}

get_average <- function(mutation_counts, column_ind = "gtexInd", tissue_pattern = "Brain", new_name = "Brain", column_average = "log2_fc_obs_expected", other_cols = c("gtexInd", "log10_stem_cell_division", "AGE")) {
    
    mutation_counts <- as.data.frame(mutation_counts[,unique(c(other_cols, column_average, "tissue"))])
    for(i in unique(mutation_counts[,column_ind])) {
        
        # the subset of data we are working with, get rid of it in original data
        this <- mutation_counts[mutation_counts[,column_ind] == i,]
        this_tissues <- this[grepl(tissue_pattern, this$tissue), ]
        
        if(nrow(this_tissues) == 0)
            next
        
        
        # get tissues we are working with
        mutation_counts <- mutation_counts[mutation_counts[,column_ind] != i,]
        this <- this[!grepl(tissue_pattern, this$tissue), ]
        
        new_append <- this_tissues[1,]
        new_append$tissue <- new_name
        new_append[,column_average] <- mean(this_tissues[,column_average])
        
        # Appeding
        this <- rbind(this, new_append)
        mutation_counts <- rbind(mutation_counts, this)
    }
    
    return(mutation_counts)
}


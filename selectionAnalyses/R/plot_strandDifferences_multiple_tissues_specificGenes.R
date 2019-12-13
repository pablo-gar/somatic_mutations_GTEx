# Plots mutation type counts by strand
# take the mutation table output from dnds script

library("ggplot2")
library("dplyr")
library("tidyr")
library("gplots")
source("../../R/ggthemes.R")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/plots.R", chdir = T)

nucleotideFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n1"
NUC_COUNTS <- read.table(nucleotideFile, sep = "\t", stringsAsFactors = F, header = F, row.names = 1)

args <- commandArgs(T)
outPrefix <- args[1] 
out_table_gwas_folder <- args[2] 
readCountFile <- args[3]  
transDiversityFile <- args[4] 
genes_list_file <- args[5]
mutationFiles <- args[-(1:5)] 

#outPrefix <- c("/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences_pairwise_comparisons/NK.cells.resting_genes/Lung/n6_0.0_0.7/")
#out_table_gwas_folder <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung/n6_0.0_0.7/"
#mutationFiles <- list.files("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout", full.names = T, recursive = T, pattern = "dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt")
#readCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"
#transDiversityFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/transcriptome_shannon_diversity.txt"
#genes_list_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/genes_Lung_0.1_NK.cells.resting.txt"

#genes_list_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/genes_Lung_0.2_NK.cells.resting.txt"
#outPrefix <- c("/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences_pairwise_comparisons/NK.cells.resting_genes/Lung_0.2/n6_0.0_0.7/")
#out_table_gwas_folder <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung_0.2/n6_0.0_0.7/"

#genes_list_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/genes_Lung_anti_0.2_NK.cells.resting.txt"
#outPrefix <- c("/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences_pairwise_comparisons/NK.cells.resting_genes/Lung_anti_0.2/n6_0.0_0.7/")
#out_table_gwas_folder <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung_anti_0.2/n6_0.0_0.7/"

#genes_list_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/genes_Lung_anti_0.1_NK.cells.resting.txt"
#outPrefix <- c("/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences_pairwise_comparisons/NK.cells.resting_genes/Lung_anti_0.1/n6_0.0_0.7/")
#out_table_gwas_folder <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/strand_fc/NK.cells.resting_genes/Lung_anti_0.1/n6_0.0_0.7/"


mutationTypes <- c(`G>A` = "C>T", `G>T` = "C>A", `G>C` = "C>G", `A>G` = "T>C", `A>C` = "T>G", `A>T` = "T>A")


#------------------------
# METHODS

read_mutation_files <- function(mutationFiles) {
    
    results <- list()
    
    for(i in seq_along(mutationFiles)) {
        current <- read.table(mutationFiles[i], sep = "\t", stringsAsFactors = F, header = T)
        current$tissue <- basename(dirname(dirname(mutationFiles[i])))
        results[[i]] <-  current
    }
    
    results <- do.call(rbind, results)
    gtexIds <- sraToGtex(results$sampleID)
    results$gtexId <- gtexIds[results$sampleID]
    
    return(results)
}

getMutationStats <- function(x, normalized = F) {
    
    mutationStats_all <- x %>%
                         #Counts mutations per person
                         group_by(gtexId, tissue, onTranscribed, mutation) %>%
                         count(mutation) 
                     
    mutationStats_all_sraId <- x %>%
                         #Counts mutations per person
                         group_by(sampleID, tissue, onTranscribed, mutation) %>%
                         count(mutation) %>%
                         ungroup() %>%
                         group_by(sampleID, tissue, mutation) %>%
                         spread(onTranscribed, n) %>%
                         mutate(fc = `TRUE`/`FALSE`) %>%
                         ungroup()
                         
                         
    if(normalized) {
        mutationStats_all$n <- mutationStats_all$n / NUC_COUNTS[ gsub("(\\w)>.+", "\\1", mutationStats_all$mutation), 1] * 2
    }
    
    mutation_paired_fc <- mutationStats_all %>%
                         #Calculates mean and stderr across all people
                         group_by(tissue, gtexId, mutation) %>%
                         spread(onTranscribed, n) %>%
                         mutate(fc = `TRUE`/`FALSE`) %>%
                         ungroup()
                         
    mutationStats <- mutationStats_all %>%
                         #Calculates mean and stderr across all people
                         group_by(onTranscribed, mutation, tissue) %>%
                         summarise(avg = mean(n), std = sd(n), stderr = sd(n) / length(n))
                     
    return(list(mutationStats,mutation_paired_fc, mutationStats_all_sraId))
                     
                         
}

    
plotScatter_and_get_cors <- function(x, outPrefix, method_p_adjust) {
    
    tissues <- unique(x$tissue)
    
    stopifnot(length(tissues) > 1)
    
    x <- x %>%
        mutate(fc = log2(fc)) %>%
        select(gtexId, tissue, mutation, fc) %>%
        spread(tissue, fc) 
    
    all_cors <- list()
    counter <- 1
    for(i in 1:(length(tissues) - 1)) {
        for(j in (i + 1) : length(tissues)) {
            #p <- scatter(as.data.frame(x), x = tissues[i], y = tissues[j], facet_x = "mutation", nrowFactor = 2, regression = T)
            #p <- p + 
            #    xlab(paste(tissues[i], "mutations log2(transcribed/non-transcribed)")) +
            #    ylab(paste(tissues[j], "mutations log2(transcribed/non-transcribed)")) +
            #    theme_noGrid()
            #ggsave(paste0(outPrefix, "_", tissues[i], "_", tissues[j], "_strand_bias_fold_change_scatter.pdf"), p)
            
            for(mut in unique(x$mutation)) {
                y <- x[x$mutation == mut,]
                
                mut_tissue_a <- y[,tissues[i], drop = T]
                mut_tissue_b <- y[,tissues[j], drop = T]
                n_inds <- sum(!is.na(mut_tissue_a) & !is.na(mut_tissue_b))
                
                if(n_inds >= 30) {
                    result_cor <- cor.test(mut_tissue_a, mut_tissue_b)
                    all_cors[[counter]] <- data.frame(mutation = mut, tissue_a = tissues[i], tissue_b = tissues[j], n_inds = n_inds, cor = result_cor[["estimate"]], pvalue = result_cor[["p.value"]], stringsAsFactors = F)
                } else {
                    all_cors[[counter]] <- data.frame(mutation = mut, tissue_a = tissues[i], tissue_b = tissues[j], n_inds = n_inds, cor = NA, pvalue = NA, stringsAsFactors = F)
                }
                    
                counter <- counter + 1
            }
        }
    }
    
    all_cors <- do.call(rbind, all_cors)
    
    # Adjusting pvalue
    all_cors <- all_cors %>% 
        group_by(mutation) %>%
        mutate(p_adjusted = p.adjust(pvalue, method = method_p_adjust)) %>%
        ungroup()
    
    return(all_cors)
    
}


make_cor_matrix <- function(x, group_column, group_value, x_column, y_column, val_column, diag_val = 0) {
    x <-  x[x[,group_column] == group_value,]
    vars <- unique(c(x[,x_column, drop = T], x[ ,y_column, drop = T]))
    
    results <- matrix(NA, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
    
    for(i in 1:nrow(x)) {
        var_a <- x[i, x_column, drop = T]
        var_b <- x[i, y_column, drop = T]
        
        results[var_a, var_b] <- results[var_b, var_a] <- x[i, val_column, drop = T]
        results[var_a, var_a] <- diag_val
        results[var_b, var_b] <- diag_val
    }
    return(results)
} 

plot_heatmap  <- function(x, group_column, group_value, x_column, y_column, val_column, diag_val = 0, pval_signif = 0.01, trunc_pval = 8) {
    
    mut_matrix <- do.call(make_cor_matrix, as.list(match.call())[2:8])
    mut_matrix[mut_matrix > trunc_pval] <- trunc_pval
    mut_matrix[mut_matrix < -trunc_pval] <- -trunc_pval
    
    cols <- get_colors_significant(trunc_pval, -log10(pval_signif))
    print(length(cols))
    to_keep <-rowSums(is.na(mut_matrix)) < (0.7 * nrow(mut_matrix))
    mut_matrix <- mut_matrix[to_keep, to_keep]
    heatmap.2(x = mut_matrix, trace = "none", density.info = "none", na.color = "grey30",
              symkey = T, col = cols, 
              breaks = seq(-trunc_pval - 0.5, trunc_pval + 0.5, length.out = length(cols) + 1)
             )  
}



write_table_for_GWAS <- function(mutationStats_all_sraId, readCountFile, transDiversityFile, out_table_gwas_folder) {
    mutationStats_all_sraId$mutation <- gsub(">", "_", mutationStats_all_sraId$mutation)
    metadata_cors <- list()
    counter <- 1
    for(currentTissue in unique(mutationStats_all_sraId$tissue)) {
        current <- mutationStats_all_sraId %>%
            filter(tissue ==  currentTissue) %>%
            select(sampleID, tissue, mutation, fc) %>%
            spread(mutation, fc) %>% 
            dplyr::rename(sample = sampleID) %>%
            select(-tissue)
        
        current <- current[,ncol(current):1]
        current <- as.data.frame(current)
        rownames(current) <- current$sample
        
        if(!grepl("EXO", currentTissue)) {
            meta <-  read_all_metadata(current$sample, readCountFile, transDiversityFile, indInfo = c("AGE", "GENDER", "BMI", "RACE"))
            current <- current[rownames(meta),]
            current <- cbind(current, meta)
        }
        write.table(current, file.path(out_table_gwas_folder, paste0(currentTissue, "-fcSrands.txt")), sep  = "\t", quote = F, col.names = T, row.names = F)
    }
    
}


#------------------------
# MAIN

# Reads genes of interest
genes_interest <- read.table(genes_list_file, sep = "\t", stringsAsFactors = F, header = T)

# Reads in mutations
mutationTable <- read_mutation_files(mutationFiles)
mutationTable <- mutationTable[mutationTable$gene %in% genes_interest[,1],]

mutationTable$mutation_cod <- paste0(mutationTable$ref_cod, ">", mutationTable$mut_cod)


# Appending if mutation is on transcribed and creating 6-type mutation profile
mutationTable$onTranscribed <- T
mutationTable[mutationTable$mutation_cod %in% mutationTypes[1:6], "onTranscribed"] <- F

mutationTable$mutation <- mutationTable$mutation_cod
mutationTable[mutationTable$onTranscribed, "mutation"] <- mutationTypes[mutationTable[mutationTable$onTranscribed, "mutation"]]


# Calculating total mutations per person
mutation_stats_nonNorm <-  getMutationStats(mutationTable, normalized = F)
cor_results <- plotScatter_and_get_cors(mutation_stats_nonNorm[[2]], outPrefix, method_p_adjust = "BH")
cor_results <- cor_results[order(cor_results$pvalue),]
cor_results$signed_p <- cor_results$cor/abs(cor_results$cor) * (-log10(cor_results$p_adjusted))

# Plots heatmaps 
#for(mut in unique(cor_results$mutation)) {
#    pdf(paste0(outPrefix, "1_strand_bias_cor_heatmap_all_tissues_", gsub(">", "_", mut), ".pdf"))
#    plot_heatmap(cor_results, "mutation", mut, "tissue_a", "tissue_b", "signed_p", diag_val = NA, pval_signif = 0.01, trunc_pval = 8)
#    dev.off()
#}

# Write correlation results
#write.table(cor_results, paste0(outPrefix, "1_cor_results.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

# Write tables for GWAS
write_table_for_GWAS (mutation_stats_nonNorm[[3]], readCountFile, transDiversityFile, out_table_gwas_folder)

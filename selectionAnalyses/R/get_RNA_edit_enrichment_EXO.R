# Plots mutation type counts by strand
# take the mutation table output from dnds script

library("ggplot2")
library("dplyr")
library("tidyr")
source("../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    #--------------------------------------------------#
    # INITIALIZING
    #--------------------------------------------------#

    #mutation_to_focus <- cmdArgs[1] #
    #original_profile <- as.logical(cmdArgs[2]) #
    #out_prefix <- cmdArgs[3]
    #mutationFile <- cmdArgs[4] #

    mutation_to_focus <- "A>G"
    original_profile <- FALSE  # Whether to work with C>T-like mutation profile or G>A-like profiles
    out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/selection/rna_edit_enrichment/Whole_Blood_EXO_expressed_genes/n6_0.0_0.7/A>G_"
    mutationFile_rna <- c("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood/n6_0.0_0.7/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt")
    mutationFile <- c("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout_EXO/Whole_Blood_EXO/n6_0.0_0.7/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt")

    mutationTypes <- c(`G>A` = "C>T", `G>T` = "C>A", `G>C` = "C>G", `A>G` = "T>C", `A>C` = "T>G", `A>T` = "T>A")
    mutationTypes_reverse <- setNames(names(mutationTypes), mutationTypes)
    
    mutation_to_focus_both <- c(mutation_to_focus, mutationTypes[mutation_to_focus], mutationTypes_reverse[mutation_to_focus])
    mutation_to_focus_both <- mutation_to_focus_both[!is.na(mutation_to_focus_both)]

    rna_edit_context <- c("AAG", "TAG", "CAG")
    
    dir.create(dirname(out_prefix), recursive = T)
    
    #--------------------------------------------------#
    # MAIN
    #--------------------------------------------------#

    # Reads in mutations
    mutationTable <- read.table(mutationFile, sep = "\t", header = T, stringsAsFactors = F)
    mutationTable_rna <- read.table(mutationFile_rna, sep = "\t", header = T, stringsAsFactors = F) 
    expressed_genes <- table(mutationTable_rna$gene)
    expressed_genes <- names(expressed_genes)
    #expressed_genes <- names(expressed_genes[expressed_genes >=3])
    
    #mutationTable <- mutationTable[paste0(mutationTable$sampleID, mutationTable$gene) %in% paste0(mutationTable_rna$sampleID, mutationTable_rna$gene),]
    mutationTable <- mutationTable[mutationTable$gene %in%  expressed_genes,]
    
    mutationTable$mutation_cod <- paste0(mutationTable$ref_cod, ">", mutationTable$mut_cod)
    
    # Keep only mutation of interest
    mutationTable <- mutationTable[ mutationTable$mutation_cod %in% mutation_to_focus_both, ]

    # Appending if mutation is on transcribed and creating 6-type mutation profile
    if(original_profile) {
        mutationTable <- get_single_mutation(mutationTable, mutationTypes)
    } else {
        mutationTable <- get_single_mutation(mutationTable, mutationTypes_reverse)
    }
    
    # Get context 
    mutationTable <- get_context(mutationTable)
    mutationTable$is_edit_context <- mutationTable$context %in% rna_edit_context
    
    
    # Plot
    counts <- mutationTable %>%
        group_by(sampleID, context) %>%
        summarise(isCoding = sum(onCoding), noCoding = sum(!onCoding), percent_coding = isCoding/(isCoding+noCoding)) %>%
        ungroup()

    counts$is_edit_context <- counts$context %in% rna_edit_context 
    unique_context <- unique(counts$context)
    counts$context <- factor(counts$context, ordered = T, 
                                     levels = c(unique_context[unique_context %in% rna_edit_context], unique_context[!unique_context %in% rna_edit_context]))
    
    
    p <- ggplot(counts, aes(x = context, y = percent_coding, fill = is_edit_context)) +
            geom_boxplot() +
            ylab("Percentage of mutations in coding strand") +
            theme_grid_y() + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

    
    # Doing total per coding in blood
    counts2 <- counts %>%
        select(context, isCoding, noCoding) %>%
        gather("coding", "count", isCoding, noCoding) %>%
        group_by(context) %>%
        mutate(total = sum(count)) %>%
        ungroup() %>%
        group_by(context, coding) %>%
        summarise(count = sum(count), percentage = count/total[1]) %>%
        ungroup()
    
        
    p_bars_count <- ggplot(counts2, aes(x = context, y = count, fill = coding)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black") +
            theme_grid_y() + 
            scale_fill_manual(values = c("#e29f0d", "#4790ef")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        
    
    p_bars_percentage <- ggplot(counts2, aes(x = context, y = percentage, fill = coding)) +
            geom_bar(stat = "identity", position = "stack") +
            ylab("Frequency") +
            theme_grid_y() + 
            scale_fill_manual(values = c("#e29f0d", "#4790ef")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        
    # Doing total 
    counts3 <- counts2 %>%
        group_by(context) %>%
        summarise(count = sum(count)) %>%
        ungroup()
    
        
    p_bars_count_all <- ggplot(counts3, aes(x = context, y = count)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "grey60") +
            theme_grid_y() + 
            scale_fill_manual(values = c("#e29f0d", "#4790ef")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
    #save plots
    ggsave(paste0(out_prefix, "coding_non_coding_boxplots.pdf"), p, height = 4)
    ggsave(paste0(out_prefix, "coding_non_coding_bar_counts.pdf"), p_bars_count, height = 4)
    ggsave(paste0(out_prefix, "coding_non_coding_bar_percentage.pdf"), p_bars_percentage, height = 4)
    ggsave(paste0(out_prefix, "all_bar_counts.pdf"), p_bars_count_all, height = 4)
        
}


#------------------------
# METHODS

get_single_mutation <- function(mutationTable, mutationTypes_this) {
                                
    mutationTable$onCoding <- F
    mutationTable[mutationTable$mutation_cod %in% mutationTypes_this, "onCoding"] <- T
    
    # Converting to single mutation
    mutationTable$mutation <- mutationTable$mutation_cod
    mutationTable[!mutationTable$onCoding, "mutation"] <- mutationTypes_this[mutationTable[!mutationTable$onCoding, "mutation"]]
    
    return(mutationTable)
    
}

get_context <- function(mutationTable){
                        
    complement <- c(`G` = "C", `C` = "G", `A` = "T", `T` = "A")

    # Get watson context
    context <- do.call(rbind, strsplit(mutationTable$ref3_cod, ""))
    context_notCod <- context[!mutationTable$onCoding,]
    context_notCod <- apply(context_notCod, 1, function(x) paste0(complement[x[3]], complement[x[2]], complement[x[1]]))

    mutationTable$context <- mutationTable$ref3_cod
    mutationTable[!mutationTable$onCoding, "context"] <- context_notCod
    
    return(mutationTable)
}

main()

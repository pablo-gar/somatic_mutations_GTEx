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

    mutation_to_focus <- cmdArgs[1] #
    original_profile <- as.logical(cmdArgs[2]) #
    out_prefix <- cmdArgs[3]
    context_counts_dir <- cmdArgs[4]
    mutationFile <- cmdArgs[5] #
    
    print(length(cmdArgs) )
    if(length(cmdArgs) >= 6) {
        mutationFile_rna <- cmdArgs[6]
    } else {
        mutationFile_rna <- NULL
    }
    

    # RNA
    #mutation_to_focus <- "A>G"
    #original_profile <- FALSE  # Whether to work with C>T-like mutation profile or G>A-like profiles
    #out_prefix <- "~/deleteme_A_G_"
    #context_counts_dir <- c("/scratch/users/paedugar/somaticMutationsProject/context_counts/1_nucleotides/Whole_Blood/")
    #mutationFile <- c("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood/n6_0.0_0.7/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt")
    #mutationFile_rna <- NULL
    
    # Colon transverse
    #mutation_to_focus <- "A>G"
    #original_profile <- FALSE  # Whether to work with C>T-like mutation profile or G>A-like profiles
    #out_prefix <- "~/deleteme_A_G_"
    #context_counts_dir <- c("/scratch/users/paedugar/somaticMutationsProject/context_counts/1_nucleotides/Colon_Transverse/")
    #mutationFile <- c("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Colon_Transverse/n6_0.0_0.7/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt")
    #mutationFile_rna <- NULL
    
    # Exome
    #mutation_to_focus <- "A>G"
    #original_profile <- FALSE  # Whether to work with C>T-like mutation profile or G>A-like profiles
    #out_prefix <- "~/deleteme_A_G"
    #context_counts_dir <- c("/scratch/users/paedugar/somaticMutationsProject/context_counts/1_nucleotides/Whole_Blood_EXO/")
    #mutationFile <- "/scratch/users/paedugar/somaticMutationsProject/selection/dndsout_EXO/Whole_Blood_EXO/n6_0.0_0.7/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt"
    #mutationFile_rna <- "/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood/n6_0.0_0.7/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt"

    mutationTypes <- c(`G>A` = "C>T", `G>T` = "C>A", `G>C` = "C>G", `A>G` = "T>C", `A>C` = "T>G", `A>T` = "T>A")
    mutationTypes_reverse <- setNames(names(mutationTypes), mutationTypes)
    nucleotides <- c("C", "T")
    nucleotides_reverse <- reverseComplement_string(nucleotides)
    
    mutation_to_focus_both <- c(mutation_to_focus, mutationTypes[mutation_to_focus], mutationTypes_reverse[mutation_to_focus])
    mutation_to_focus_both <- mutation_to_focus_both[!is.na(mutation_to_focus_both)]
    
    if (original_profile) {
        mutation_base <- "T"
    } else {
        mutation_base <- "A"
    }
        

    rna_edit_context <- c("AAG", "TAG", "CAG")
    
    #--------------------------------------------------#
    # MAIN
    #--------------------------------------------------#

    #--------------------------------------------------#
    # Reads contexts
    context_counts <- read_context_counts(context_counts_dir)
    context_counts <- context_counts[!is.na(context_counts$total),]
    if(original_profile) {
        context_counts <- get_coding(context_counts, nucleotides)
    } else {
        context_counts <- get_coding(context_counts, nucleotides_reverse)
    }
    context_counts_global <- context_counts %>%
        group_by(context_final, coding) %>%
        summarise(count = sum(total)) %>%
        ungroup()
    
    context_counts_global$is_edit_context <- context_counts_global$context_final %in% rna_edit_context 
    
    # Keep only mutation of interest
    context_length <- nchar(context_counts_global$context_final[1])
    middle_base <- ceiling(context_length/2)
    middle_nucs <- sapply(strsplit(context_counts_global$context_final, ""), function(x) x[middle_base])
    context_counts_global <- context_counts_global[ middle_nucs %in% mutation_base, ]
    
    # converting to ordered factors
    unique_context <- unique(context_counts_global$context_final)
    context_counts_global$context_final <- factor(context_counts_global$context_final, ordered = T, 
                                     levels = c(unique_context[unique_context %in% rna_edit_context], unique_context[!unique_context %in% rna_edit_context]))
    
    #--------------------------------------------------#
    # Reads in mutations
    mutationTable <- read.table(mutationFile, sep = "\t", header = T, stringsAsFactors = F)
    
    
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
    
    # In case we are working with exome file
    if(!is.null(mutationFile_rna)) {
        mutationTable_rna <- read.table(mutationFile_rna, sep = "\t", header = T, stringsAsFactors = F) 
        expressed_genes <- table(mutationTable_rna$gene)
        expressed_genes <- names(expressed_genes)
        
        mutationTable <- mutationTable[mutationTable$gene %in%  expressed_genes,]
        mutationTable$is_edit_context <- mutationTable$context %in% rna_edit_context
    }
    
        
    
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
        ungroup() %>%
        as.data.frame()
    
    # append background
    context_counts_global <- as.data.frame(context_counts_global)
    rownames(context_counts_global) <- paste0(context_counts_global$context, ".", context_counts_global$coding)
    
    counts2$background_count <- context_counts_global[paste0(counts2$context, ".", counts2$coding), "count"]
    counts2$count_norm <- counts2$count/ counts2$background_count
        
    p_bars_count <- ggplot(counts2, aes(x = context, y = count, fill = coding)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black") +
            theme_grid_y() + 
            scale_fill_manual(values = c("#e29f0d", "#4790ef")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
    p_bars_count_norm <- ggplot(counts2, aes(x = context, y = count_norm, fill = coding)) +
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
        
    p_bars_count_background <- ggplot(context_counts_global, aes(x = context_final, y = count, fill = coding)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black") +
            theme_grid_y() + 
            scale_fill_manual(values = c("#e29f0d", "#4790ef")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        
    # Doing total 
    counts3 <- counts2 %>%
        group_by(context) %>%
        summarise(count = sum(count)) %>%
        ungroup()
    
    context_counts_global_total <- context_counts_global %>%
        group_by(context_final) %>%
        summarise(count = sum(count)) %>%
        ungroup()
    
        
    p_bars_count_all <- ggplot(counts3, aes(x = context, y = count)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "grey60") +
            theme_grid_y() + 
            scale_fill_manual(values = c("#e29f0d", "#4790ef")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
    p_bars_count_all_background <- ggplot(context_counts_global_total, aes(x = context_final, y = count)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "grey60") +
            theme_grid_y() + 
            scale_fill_manual(values = c("#e29f0d", "#4790ef")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
    #save plots
    ggsave(paste0(out_prefix, "coding_non_coding_boxplots.pdf"), p, height = 4, width = 9)
    ggsave(paste0(out_prefix, "coding_non_coding_bar_counts.pdf"), p_bars_count, height = 4, width = 9)
    ggsave(paste0(out_prefix, "coding_non_coding_bar_counts_background.pdf"), p_bars_count_background, height = 4, width = 9)
    ggsave(paste0(out_prefix, "coding_non_coding_bar_counts_normalized.pdf"), p_bars_count_norm, height = 4, width = 9)
    ggsave(paste0(out_prefix, "coding_non_coding_bar_percentage.pdf"), p_bars_percentage, height = 4, width = 9)
    
    ggsave(paste0(out_prefix, "all_bar_counts.pdf"), p_bars_count_all, height = 4, width = 9)
    ggsave(paste0(out_prefix, "all_bar_counts_background.pdf"), p_bars_count_all_background, height = 4, width = 9)
        
}


#------------------------
# METHODS

read_context_counts <- function(context_counts_dir) {
    
    results <- list()
    
    strand_dirs <- list.dirs(context_counts_dir)
    stopifnot(any(grepl("strand_-", strand_dirs)), any(grepl("strand_\\+", strand_dirs)))
    
    
    # Reading data
    i <- 1
    for(strand_dir in c("strand_-", "strand_+")) {
        count_files <- list.files(file.path(context_counts_dir, strand_dir), full = T)
        for(count_file in count_files) {
            current <- read.table(count_file, sep = "\t", header = F, stringsAsFactors = F)
            colnames(current) <- c("context", "counts")
            current$sample <- gsub(".txt", "", basename(count_file))
            current$strand <- strand_dir
            results[[i]] <- current
            i <- i + 1
        }
    }
    
    # Gathering data
    results <- do.call(rbind, results)
    results$strand <- gsub("strand_-", "minus", results$strand)
    results$strand <- gsub("strand_\\+", "plus", results$strand)
    
    # Reverse complememtns context in negative strand
    results[results$strand == "minus", "context"] <- reverseComplement_string(results[results$strand == "minus", "context"])
    
    # Re format to wide 
    results <- spread(results, strand, counts)
    results$total <- results$minus + results$plus
    
    
    return(results)
   
    
}

reverseComplement_string <- function(x, Reverse = T) {
    
    x <- toupper(x)
    key <- c(A = "T", T = "A", C = "G", G = "C", N = "N")
    x <- strsplit(x, "")
    
    result <- sapply(x, function(x) { if(Reverse) x <- rev(x); paste(key[x], collapse = "") })
    
    return(result)
    
}

get_single_mutation <- function(mutationTable, mutationTypes_this) {
                                
    mutationTable$onCoding <- F
    mutationTable[mutationTable$mutation_cod %in% mutationTypes_this, "onCoding"] <- T
    
    # Converting to single mutation
    mutationTable$mutation <- mutationTable$mutation_cod
    mutationTable[!mutationTable$onCoding, "mutation"] <- mutationTypes_this[mutationTable[!mutationTable$onCoding, "mutation"]]
    
    return(mutationTable)
    
}

get_coding <- function(context_counts, nucleotide_this = c("A", "G")) {
                                
    context_counts$coding <- F
    context_length <- nchar(context_counts$context[1])
    middle_base <- ceiling(context_length/2)
    
    nucleotides <- sapply(strsplit(context_counts$context, ""), function(x) x[middle_base])
    
    context_counts[nucleotides %in% nucleotide_this, "coding"] <- T
    
   # Getting reverse complements of non-coding
    context_counts$context_final <- context_counts$context
    context_counts[!context_counts$coding, "context_final"] <- reverseComplement_string(context_counts[!context_counts$coding, "context_final"])
    
    context_counts$coding <- ifelse(context_counts$coding, "isCoding", "noCoding")
    
    return(context_counts)
    
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

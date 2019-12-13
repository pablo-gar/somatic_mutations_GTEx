# Takes mutation files with two experiments and counts how many times mutations
# are found in both experiments

# Usage
# Rscript get_agreement.R /scratch/users/paedugar/somaticMutationsProject/method_validation/mutations/RNA_first outputTable.txt

library("dplyr")
library("ggplot2")
source("../../R/ggthemes.R")
source("../../R/misc.R", chdir = T)

##------------------------------------------------------------
## GLOBAL
nAtLeastCovered <- 8

# Indeces for reference and alternate alleles counts, for each experiment 
REF_I <- c(7, 12) 
ALT_I <- c(6, 11) 
COV_I <- c(5, 10) 

# Indeces for reference and alternate allele names
REF_A <- c(3, 8)
ALT_A <- c(4, 9)

# Indeces for reference context and strand
CTXT_I <- 13
STRAND_I <- 14

# Colors for false-positves and true-positives
boxColors <- c("#f2aa82",  "#448cff")

rna_edit_context <- c("AAG", "TAG", "CAG")

##------------------------------------------------------------

main <- function(cmdArgs = commandArgs(T)) {
    
    mutationDir <- cmdArgs[1]
    outPrefix <- cmdArgs[2]

    #mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/method_validation/mutations_lax/RNA_first_context/" 
    #outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/method_validation/fdr/RNA_first_context/rna_edit_"
    dir.create(dirname(outPrefix), recursive = T)
    
    # Read mutations
    tableFile <- getTableFile(mutationDir)
	mutations <- getMutations(tableFile)
    mutations <- do.call(rbind, mutations)
    
    # Select rna edit-like mutations
    mutations <- mutations[ (mutations[,REF_A[1]] == "T" & mutations[,ALT_A[1]] == "C") | (mutations[,REF_A[1]] == "A" & mutations[,ALT_A[1]] == "G") ,]
    
    mutations <- get_coding(mutations)
    #mutations <- get_context_mut(mutations)
    
    mutations <- removeLowCovered(mutations, n = nAtLeastCovered)
    
    # Getting FDRs per context and coding
    fdr <- mutations %>%
        group_by(V13, coding) %>%
        summarise(total = n(), realPos = sum(V4 == V9), falsePos = total - realPos, fdr = falsePos/total) %>%
        ungroup()
    
    unique_context <- unique(fdr$V13)
    fdr$context <- factor(fdr$V13, ordered = T, 
                                     levels = c(unique_context[unique_context %in% rna_edit_context], unique_context[!unique_context %in% rna_edit_context]))
    
    # Getting average fdr per sample
    fdr_average <- mutations %>%
        group_by(V13, coding, ind) %>%
        summarise(total = n(), realPos = sum(V4 == V9), falsePos = total - realPos, fdr = falsePos/total) %>%
        ungroup() %>% 
        group_by(V13, coding) %>%
        summarise(total_mean = mean(total), fdr_mean = mean(fdr),
                  total_low = bootstrap_confidence_interval(total, mean, 10e3)[1], total_high = bootstrap_confidence_interval(total, mean, 10e3)[2],
                  fdr_low = bootstrap_confidence_interval(fdr, mean, 10e3)[1], fdr_high = bootstrap_confidence_interval(fdr, mean, 10e3)[2]) %>%
        ungroup()
    
    unique_context <- unique(fdr_average$V13)
    fdr_average$context <- factor(fdr_average$V13, ordered = T, 
                                     levels = c(unique_context[unique_context %in% rna_edit_context], unique_context[!unique_context %in% rna_edit_context]))
    
    # Plotting
    p_fdr_all <- ggplot(fdr, aes(x = context, y = fdr * 100)) +
        geom_bar(aes(fill = coding), stat = "identity", position = "dodge", colour = "black") +
        xlab("Context") +
        ylab("FDR") +
        theme_grid_y()
    
    p_counts_all <- ggplot(fdr, aes(x = context, y = total)) +
        geom_bar(aes(fill = coding), stat = "identity", position = "dodge", colour = "black") +
        xlab("Context") +
        ylab("Total mutations") +
        theme_grid_y()
        
    p_fdr_average <- ggplot(fdr_average, aes(x = context, y = fdr_mean * 100, fill = coding)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
        geom_errorbar(aes(ymin = fdr_low * 100, ymax = fdr_high * 100), position = position_dodge(width = 0.9), width = 0.4) +
        xlab("Context") +
        ylab("Average FDR") +
        theme_grid_y()
    
    p_counts_average <- ggplot(fdr_average, aes(x = context, y = total_mean, fill = coding)) +
        geom_bar(stat = "identity", position = "dodge", colour = "black") +
        geom_errorbar(aes(ymin = total_low, ymax = total_high), position = position_dodge(width = 0.9), width = 0.4) +
        xlab("Context") +
        ylab("Average mutations") +
        theme_grid_y()
        
    ggsave(paste0(outPrefix, "fdr_all.pdf"), p_fdr_all, height = 4, width = 9)
    ggsave(paste0(outPrefix, "counts_all.pdf"), p_counts_all, height = 4, width = 9)
    ggsave(paste0(outPrefix, "fdr_sample_average.pdf"), p_fdr_average, height = 4, width = 9)
    ggsave(paste0(outPrefix, "counts_sample_average.pdf"), p_counts_average, height = 4, width = 9)
    
    
}

reverseComplement_string <- function(x, Reverse = T) {
    
    x <- toupper(x)
    key <- c(A = "T", T = "A", C = "G", G = "C", N = "N", X = "X")
    x <- strsplit(x, "")
    
    result <- sapply(x, function(x) { if(Reverse) x <- rev(x); paste(key[x], collapse = "") })
    
    return(result)
    
}

#' Get mutation switches in coding sequences
get_coding <- function(x, nucleotide_this = c("A", "G")) {
    
    # First convert mutations to coding side
    for(j in c(REF_A, ALT_A, CTXT_I)) 
        x[x[,STRAND_I] == "-", j] <- reverseComplement_string(x[x[,STRAND_I] == "-", j])
    
    
    # Simplify mutation code
    x$coding <- F
    x[x[,REF_A[1]] %in% nucleotide_this, "coding"] <- T
    
   # Getting reverse complements of non-coding
    for(j in c(REF_A, ALT_A, CTXT_I)) 
        x[!x$coding, j] <- reverseComplement_string(x[!x$coding, j])
    
    x$coding <- ifelse(x$coding, "isCoding", "noCoding")
    
    return(x)
    
    # Then compress from T>C/A>G to A>G
    #x$coding <- F
    #context_length <- nchar(x[1,CTXT_I])
    #middle_base <- ceiling(context_length/2)
    
}

get_context_mut <- function(x) {
    
    context_length <- nchar(x[1,CTXT_I])
    middle_base <- ceiling(context_length/2)
    contexts <- do.call(rbind, strsplit(x[,CTXT_I], ""))
    
    for(j in c(REF_A, ALT_A)) {
        current_contexts <- contexts
        current_contexts[,middle_base] <- x[, j]
        x[, j] <- apply(current_contexts, 1, paste, collapse = "")
    }
    
    return(x)
}


getTableFile <- function(x) {
    
    # Returns a data frame where each columns has info for 
    # each pairwsie-pileup comparison file found in the pileupDir
    #
    # @args x - string - path to location of pileup files
    #
    # @return - data.frame - containing the info of each pileup file:
    #     - individual id
    #     - exp A id
    #     - exp B id
    #     - coverage cutoff
    #     - nSupporting reads cutoff
    #     - path to file
    
    pileupFiles <- list.files(x)
    pileupFiles <- pileupFiles[grep(".txt$", pileupFiles)]
    
    results <- data.frame(ind = "", sampleA = "", sampleB = "", coverage = 0, nSupporting = 0, location = file.path(x, pileupFiles), stringsAsFactors = F)
                          
    pileupFiles <- gsub(".txt", "", pileupFiles)
    pileupFiles <- strsplit(pileupFiles, "-")
    
    for (i in 1:length(pileupFiles)) {
        mut_settings <- unlist(strsplit(pileupFiles[[i]][1], "_"))
        currentCoverage <-  paste(mut_settings[2:3], collapse = "_")
        currentNsupporting <-  as.integer(gsub("n", "", mut_settings[1]))
        
        results[i, "ind"] <- pileupFiles[[i]][2]
        results[i, "sampleA"] <- pileupFiles[[i]][2]
        results[i, "sampleB"] <- pileupFiles[[i]][3]
        results[i, "coverage"] <- currentCoverage
        results[i, "nSupporting"] <- currentNsupporting
    }
    
    return(results)
}


getMutations <- function(tableFile) {
    
    # Reads all pileups from the tableFile into a list of data.frames
    # appends to the pileups the coverage and nSupporting cutoffs as well as
    # individual id and experiments ids

    allPileups <- list()
    i <- 1
    for(coverage in unique(tableFile$coverage)) {
        for(nSupporting in unique(tableFile$nSupporting)){
            
            currentFiles <- tableFile[ tableFile$coverage == coverage & tableFile$nSupporting == nSupporting,]
            
            for(j in 1:nrow(currentFiles)) {
                
                currentFile <- currentFiles[j, "location"]
                currentPileup <- read.table(currentFile, sep = "\t", header = F, stringsAsFactors = F)

                # Appending info
                currentPileup$coverage <- coverage
                currentPileup$nSupporting <- nSupporting
                currentPileup$ind <- currentFiles[j, "ind"]
                currentPileup$sampleA <- currentFiles[j, "sampleA"]
                currentPileup$sampleB <- currentFiles[j, "sampleB"]

                allPileups[[i]] <- currentPileup

                i <- i + 1
            }
        }
    }

    return(allPileups)
}


removeLowCovered <- function(x, n = 5) {
    
    # From 2-experiment pileup table remove positions for which in the second experiment (Exome) we don't enough coverage
    # to observe the variant from RNA. For example a variant from RNA that was covered 2000x
    # and only observed 4x, means that we would expected 4/2000 in exome. If the same
    # position was covered 40 times in exome, then we expect (4/2000)*40 reads having the variant
    # that is 0.08 reads
    #
    # @param x - data.frame - pileup  from  readPileup functions, 11 columns
    # @param n - int - lower bound threshold for number of reads that we expect to observe in second experiment
    #
    # @return - data.frame - same pileup table as input without the filtered positions
    
    expected <- (x[,ALT_I[1]] / x[,COV_I[1]] ) * x[,COV_I[2]]
    x <- x[ expected >= n, ]
    
    return(x)
}

main()

# Takes mutation files with two experiments and counts how many times mutations
# are found in both experiments

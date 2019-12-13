# Takes mutation files with two experiments and counts how many times mutations
# are found in both experiments

# Usage
# Rscript get_agreement.R /scratch/users/paedugar/somaticMutationsProject/method_validation/mutations/RNA_first outputTable.txt

library("dplyr")
library("tidyr")
library("ggplot2")
source("../../R/ggthemes.R")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/plots.R", chdir = T)

##------------------------------------------------------------
## GLOBAL
nAtLeastCovered <- 8

MUT_KEY <- matrix(ncol = 2, byrow = T,
                  data = c(
                    "C.T", "G.A",
                    "C.A", "G.T",
                    "C.G", "G.C",
                    "T.A", "A.T",
                    "T.C", "A.G",
                    "T.G", "A.C"
                  )
                )

# Indeces for reference and alternate alleles counts, for each experiment 
REF_I <- c(7, 12) 
ALT_I <- c(6, 11) 
COV_I <- c(5, 10) 

# Indeces for reference and alternate allele names
REF_A <- c(3, 8)
ALT_A <- c(4, 9)

# Colors for false-positves and true-positives
boxColors <- c("#f2aa82",  "#448cff")

##------------------------------------------------------------

main <- function(cmdArgs = commandArgs(T)) {
    
    mutationDir <- cmdArgs[1]
    outPrefix <- cmdArgs[2]

    #mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/method_validation/mutations/RNA_first/" 
    #mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/method_validation_with_repeated/mutations_lax/RNA_first/" 
    #mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/method_validation/mutations_lax/RNA_first/" 
    #readCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"
    #transDiversityFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/transcriptome_shannon_diversity.txt"
    #outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/method_validation/fdr/RNA_first/fdr"
    #freq_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/all_mutations_frequency/mutation_freq_table_withRepeated.txt"
    
    tableFile <- getTableFile(mutationDir)
	mutations <- getMutations(tableFile)
    
    ### TO INCORPORATE mutation frequencies from all tissues
    #
    #freq_table <- read.table(freq_table_file, sep = "\t", stringsAsFactors = F, header = F)
    #mutations <- append_mut_freq(mutations, freq_table)
    #
    #
    ## Calculating maf and frequency in all tissues
    ##freq <- unlist(lapply(mutations, function(x) x$freq))
    ##maf <- unlist(lapply(mutations, function(x) x$V6 / x$V5))
    ##
    ##median_freq <- median(freq)
    ##median_maf <- median(maf)
    ##
    ##mutations <- lapply(mutations, function(x) x[x$freq > median_freq,])
    ##mutations <- lapply(mutations, function(x) x[x$freq >= 0.2,])
    
    #mutations <- lapply(mutations, function(x) x[x$freq <= 0.03,])
    
    #mutations <- lapply(mutations, function(x) x[x$V6/ x$V5 > median_maf,])
    
    #for(permute in c(F,T)) {
    for(permute in c(F,T)) {
        countsResults <- getCountTable(mutations, permute)
        alt_counts_freq <- countsResults[[2]]
        maf_freq <- countsResults[[3]]
        
        
        fdr <- get_FDR_byMutbySample(alt_counts_freq)
        fdr <- fdr[fdr$totalByMut != 0,]
        
        # Get table for histogram of n reads supporting variant alelleA
        freq_table <- alt_counts_freq %>%
            group_by(Var1) %>%
            summarise(total = sum(total))
        freq_table$total <-  100 * freq_table$total / sum(freq_table$total)
        p_freq <- ggplot(freq_table, aes(x = Var1, y = total)) +
            geom_bar(stat = "identity") +
            xlab("Reads supporting alternate") +
            ylab("Percentage") +
            coord_cartesian(xlim = c(0, 40), ylim = c(0,100)) +
            theme_grid_y()
        
        
        #Getting complement
        fdr$mutation <- setNames(c("total", "T.G", "T.C", "T.A", "C.A", "C.G", "C.T", "C.T", "C.G", "C.A", "T.A", "T.C", "T.G"),
                                 c("total", "T.G", "T.C", "T.A", "G.T", "G.C", "G.A", "C.T", "C.G", "C.A", "A.T", "A.G", "A.C"))[fdr$mutation]
        
        fdr_total <- get_FDR_total(alt_counts_freq)
        fdr_median <- get_FDR_conf_interval(fdr, mean)
        
        countsPerMut <- get_total_countsByMut(fdr)
        
        p_byMut  <- plotFDR_byMut(fdr, countsPerMut)
        p_all <- plotFDR_all(fdr)
        
        y_max <- ifelse(max(fdr_median$percentage_mean) >= 45, 100, 45)
        p_median <- ggplot(fdr_median, aes(x = mutation, y = percentage_mean)) +
            geom_bar(stat = "identity", width = 0.8, fill = rev(c("grey50", "#F1CAC9",  "#9CD169",  "#CAD0CE", "#FC3218","#2B2E34","#3BC9F3"))) + 
            geom_errorbar(aes(ymax = up_conf, ymin = low_conf), width = 0.15) +
            ylab("Mean FDR") +
            ylim(0,y_max) +
            coord_flip() +
            theme_grid_x()

        ggsave(paste0(outPrefix, "_byMut_permute_", permute, ".pdf"), p_byMut)
        ggsave(paste0(outPrefix, "_all_permute_", permute, ".pdf"), p_all, width = 4)
        ggsave(paste0(outPrefix, "_fdr_median_conf_interval_permute_", permute, ".pdf"), p_median, height = 5)
        ggsave(paste0(outPrefix, "_hist_n_reads_support_alt", permute, ".pdf"), p_freq, width = 7, height = 5)
        write.table(fdr_total, paste0(outPrefix, "_1_fdr_", permute, ".txt"), sep = "\t", quote = F, row.names = F, col.names = T)
        write.table(fdr_median, paste0(outPrefix, "_1_fdr_median", permute, ".txt"), sep = "\t", quote = F, row.names = F, col.names = T)
        
        
        # UNCOMENT THIS IF YOU WANT TO CORRELETA WITH METADATA
        #fdr_allMuts <- fdr[ fdr$mutation == "total" & fdr$reads_support == "false_positives", ]
        #metadata <- read_all_metadata(fdr_allMuts$ind, readCountFile = readCountFile, transDiversityFile = transDiversityFile, metadataCols =  c("SMRIN", "SMTSISCH", "SMTSPAX"), indInfo = c("AGE", "GENDER", "BMI"))
        #metadata <- metadata[rownames(metadata) %in% fdr_allMuts$ind,]
        #fdr_allMuts <- fdr_allMuts[fdr_allMuts$ind %in% rownames(metadata),]
        #metadata$fdr <- 0
        #metadata[fdr_allMuts$ind, "fdr"] <- fdr_allMuts$percentage
    }
}

#' Append the mutation frequency across all tissues
append_mut_freq <- function(mutations, freq_table) {
    
    rownames(freq_table) <- paste0(freq_table$V4, ".", freq_table$V5, ".", freq_table$V6, ".", freq_table$V7)
    
    for(i in seq_along(mutations)) {
        
        mutations[[i]]$freq <- NA
        ids <- paste0(mutations[[i]][,1], ".", mutations[[i]][,2], ".", mutations[[i]][,3], ".", mutations[[i]][,4])
        mutations[[i]]$freq <- freq_table[ids,1]
    }
    
    return(mutations)
}

get_FDR_byMutbySample <- function(alt_counts_freq) {
    
    alt_counts_freq %>%
    mutate(reads_support = ifelse(Var1 == 0, "false_positives", "true_positives")) %>%
    select(-Var1) %>%
    group_by(ind, reads_support) %>%
    summarise_at(1:13, sum) %>%
    ungroup() %>%
    gather("mutation", "count", C.T:total) %>%
    spread(reads_support, count) %>%
    mutate(totalByMut = true_positives + false_positives) %>%
    gather("reads_support", "count", false_positives:true_positives) %>%
    mutate(percentage = count / totalByMut * 100)
    
}

get_FDR_total <- function(alt_counts_freq) {
    
    alt_counts_freq %>%
    mutate(reads_support = ifelse(Var1 == 0, "false_positives", "true_positives")) %>%
    select(-Var1) %>%
    group_by(reads_support) %>%
    summarise_at(1:13, sum) %>%
    ungroup() %>%
    gather("mutation", "count", C.T:total) %>%
    spread(reads_support, count) %>%
    mutate(totalByMut = true_positives + false_positives) %>%
    gather("reads_support", "count", false_positives:true_positives) %>%
    mutate(percentage = count / totalByMut * 100)
    
}

get_FDR_conf_interval <- function(fdr, FUN) {
    
    fdr %>% 
        filter(reads_support == "false_positives", !is.na(percentage)) %>% 
        group_by(mutation) %>% 
        summarise(percentage_mean = FUN(percentage), low_conf = bootstrap_confidence_interval(percentage, FUN, 10e3)[1], up_conf = bootstrap_confidence_interval(percentage, FUN, 10e3)[2]) %>% 
        ungroup()
    
}


get_total_countsByMut <- function (alt_counts_freq_total) {
    alt_counts_freq_total %>%
    filter(reads_support == "false_positives") %>% 
    group_by(mutation) %>% 
    summarise(total = paste0("n=", sum(totalByMut)), reads_support =  "false_positives", percentage = Inf) %>% 
    ungroup()
}

plotFDR_byMut <- function(alt_counts_freq_total, countsPerMut) {

    ggplot(filter(alt_counts_freq_total, mutation != "total") , aes(x = reads_support, y = percentage, fill = reads_support)) +
        geom_violin() +
        geom_boxplot(width = 0.3) +
        geom_text(aes(label = total), data = filter(countsPerMut, mutation != "total"), vjust = 1) +
        facet_wrap(~mutation, ncol = 4) +
        scale_fill_manual(values = boxColors) +
        ylab("Percentage") + xlab("") + 
        theme_grid_y() +
        theme(legend.position = "top", legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plotFDR_all <- function(alt_counts_freq_total) {
    ggplot(filter(alt_counts_freq_total, mutation == "total") , aes(x = reads_support, y = percentage, fill = reads_support)) +
        geom_violin() +
        geom_boxplot(width = 0.3) +
        annotate("text", 
                 y = max(filter(alt_counts_freq_total, mutation == "total" & reads_support == "false_positives")$percentage + 2), 
                 x = "false_positives", 
                 label = paste("Mean = ", round(mean(filter(alt_counts_freq_total, mutation == "total" & reads_support == "false_positives")$percentage), 2)), 
                 vjust = 0, colour = boxColors[1], fontface = "bold") +
        annotate("text", 
                 y = min(filter(alt_counts_freq_total, mutation == "total" & reads_support == "true_positives")$percentage - 2), 
                 x = "true_positives", 
                 label = paste("Mean = ", round(mean(filter(alt_counts_freq_total, mutation == "total" & reads_support == "true_positives")$percentage), 2)), 
                 vjust = 1, colour = boxColors[2], fontface = "bold") +
        scale_fill_manual(values = boxColors) +
        ylab("Percentage") + xlab("") + 
        theme_grid_y() +
        theme(legend.position = "top", legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
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

getCountTable <- function(pileupList, permute = F) {

    # Takes a list of pileup data.frames (output of getPileups()) and returns
    # a SNV count table with the following info
    #
    # @return - data.frame - 7 columns:
    #    - Total of SNVs occurring in experiment 1 at nSupproting cutoff (all calls)
    #    - Total of SNVs occurring in both experiments at nSupproting cutoff (confident calls)
    #    - SNVs agreeing in both experiments regardless of directionality (e.g A>C and C>A) in confident calls
    #    - SNVs agreeing in both experiments in directionality (e.g A>C and A>C) in confident calls
    #    - SNVs agreeing in both experiments regardless of directionality (e.g A>C and C>A) in all calls
    #    - SNVs agreeing in both experiments in directionality (e.g A>C and A>C) in all calls
    #    - MUTATION TYPE


    allCounts <- list()
    alt_counts_freq <- list()
    maf_freq <- list()

    for (i in 1:length(pileupList)) {

        currentPileup <- pileupList[[i]]

        # Getting mutation counts
        results <- getCountsByType(currentPileup, nSupporting, permute = permute)
        currentCount <- results[[1]]
        alt_counts_freq[[i]] <- results[[2]]
        maf_freq[[i]] <- results[[3]]
        
        # Merge analogous mutations
        currentCount <- simplifyMutationTypes(currentCount)

        # Appending info
        currentCount$coverage <- currentPileup$coverage[1]
        currentCount$nSupporting <- currentPileup$nSupporting[1]
        currentCount$ind <- currentPileup$ind[1]
        currentCount$sampleA <- currentPileup$sampleA[1]
        currentCount$sampleB <- currentPileup$sampleB[1]
        
        alt_counts_freq[[i]]$ind <- currentPileup$ind[1]
        maf_freq[[i]]$ind <- currentPileup$ind[1]


        allCounts[[i]] <-  currentCount
    }

    
    allCounts <- do.call(rbind, allCounts)
    alt_counts_freq <- do.call(rbind, alt_counts_freq)
    maf_freq <- do.call(rbind, maf_freq)

    return(list(allCounts, alt_counts_freq, maf_freq))
}

getCountsByType <- function(x, nSupporting, permute = F) {

    # Gets different count values from a pileup file with 2 experiments
    # Gets counts separated by MUTATION TYPE
    #
    # @args x - data.frame - pileup  from  readPileup functions, 11 columns
    # @args nSupporting - int - cutoff: number of read supporting the variant
    #
    # @return - data.frame - 7 columns:
    #    - Total of SNVs occurring in experiment 1 at nSupproting cutoff (all calls)
    #    - Total of SNVs occurring in both experiments at nSupproting cutoff (confident calls)
    #    - SNVs agreeing in both experiments regardless of directionality (e.g A>C and C>A) in confident calls
    #    - SNVs agreeing in both experiments in directionality (e.g A>C and A>C) in confident calls
    #    - SNVs agreeing in both experiments regardless of directionality (e.g A>C and C>A) in all calls
    #    - SNVs agreeing in both experiments in directionality (e.g A>C and A>C) in all calls
    #    - MUTATION TYPE


    results <- list()
    alt_count_freq <- maf_freq <- lapply(as.vector(MUT_KEY), function(x) as.table(c(`0`= 0)))
    names(alt_count_freq) <- names(maf_freq) <- as.vector(MUT_KEY)
    
    mutTypes <- paste0(x[, REF_A[1]], ".", x[, ALT_A[1]])

    for (mutType in unique(mutTypes)) {
        currentTable <- x[mutTypes == mutType, ]
        nSupporting <- currentTable$nSupporting[1]
        
        countResults <- getCounts(currentTable, nSupporting, permute = permute)
        
        # Summary table with counts
        currentCounts <- countResults[[1]]
        # Frequency tables of alt allele counts and MAFs
        alt_count_freq[[mutType]]  <- countResults[[2]]
        maf_freq[[mutType]]  <- countResults[[3]]

        results[[mutType]] <- currentCounts
    }
    

    # Merging results
    results <- as.data.frame(do.call(rbind, results))
    results$mutType <- unique(mutTypes)
    rownames(results) <- 1:nrow(results)
    
    
    # Merging frequency tables into data.frames
    alt_count_freq2 <- alt_count_freq
    alt_count_freq <- mergeFrequencyTables(alt_count_freq2)
    maf_freq <- mergeFrequencyTables(maf_freq)
    


    return(list(results, alt_count_freq, maf_freq))
}

mergeFrequencyTables <- function(x, isNumeric = T) {
    
    # x is a list of outputs from table()
    # returns a data.frame where each col is an element of x
    
    col.names = names(x)
    
    merge1 <- function(x, y) merge(x, y, by = 1, all = TRUE)
    out <- Reduce(merge1, lapply(x, as.data.frame, stringsAsFactors = FALSE))
    colnames(out)[-1] <- col.names
    out[is.na(out)] <- 0
    out$total <- rowSums(out[,-1])
    if(isNumeric) {
        out[,1] <- as.numeric(out[,1])
        out <- out[order(out[,1]),]
    }
    return(out)
}
    



getCounts <- function(x, nSupporting, permute = F) {
   
    # Gets different count values from a pileup file with 2 experiments
    # 
    # @args x - data.frame - pileup  from  readPileup functions, 11 columns
    # @args nSupporting - int - cutoff: number of read supporting the variant
    #
    # @return - data.frame - 7 columns:
    #    - Total of SNVs occurring in experiment 1 at nSupproting cutoff (all calls)
    #    - Total of SNVs occurring in both experiments at nSupproting cutoff (confident calls)
    #    - SNVs agreeing in both experiments regardless of directionality (e.g A>C and C>A) in confident calls
    #    - SNVs agreeing in both experiments in directionality (e.g A>C and A>C) in confident calls
    #    - SNVs agreeing in both experiments regardless of directionality (e.g A>C and C>A) in all calls
    #    - SNVs agreeing in both experiments in directionality (e.g A>C and A>C) in all calls
   
    countTypes <- c("Total_SNV.confident",
                   "SNV_agreed_bases.confident",
                   "SNV_agreed_direction.confident",
                   "SNV_disagreed.confident",
                   "Total_SNV.all",
                   "SNV_agreed_bases.all",
                   "SNV_agreed_direction.all",
                   "SNV_disagreed.all"
                   )
   
    counts <- rep(0, length(countTypes))
    names(counts) <- countTypes
   
    # Not count positions for which in the second experiment (Exome) we don't enough coverage
    # to observe the variant from RNA. For example a variant from RNA that was covered 2000x
    # and only observed 4x, means that we would expected 4/2000 in exome. If the same
    # position was covered 40 times in exome, then we expect (4/2000)*40 reads having the variant
    # that is 0.08 reads
    x <- removeLowCovered(x, n = nAtLeastCovered)
    
    # If no mutations return base objects
    if(nrow(x) == 0) {
        counts <- counts
        alt_count_freq <- as.table(c(`0` = 0))
        maf_freq <- as.table(c(`0` = 0))
    } else {
        
        # Permute mutations if desired 
        if(permute) {
            alt <- x[1,ALT_A[1]]
            ref <- x[1,REF_A[1]]
            
            bases <- c("A", "T", "G", "C")
            bases <- bases[ !bases %in% c(alt, ref)  ]
            #bases <- bases[ !bases %in% c(ref)  ]
            
            x[,ALT_A[1]] <- sample(bases, nrow(x), replace = T)
            #x[x[,ALT_A[2]] != "X",ALT_A[2]] <- sample(bases, sum(x[,ALT_A[2]] != "X"), replace = T)
        }
       
        # Merging ref and alt allele names for easier counting
        confident <- x[, ALT_I[2]] >= nSupporting
        mutA <- paste0(x[,REF_A[1]], x[,ALT_A[1]])
        mutB <- paste0(x[,REF_A[2]], x[,ALT_A[2]])
        noCallMutB <- grepl("[Xx]", mutB) # Only one call in second experiment
       
        # Invert names in one experiment to look for inversion calls
        mutBinverted <- paste0(x[,ALT_A[2]], x[,REF_A[2]])
       
        # Counting
        # Total muts
        counts["Total_SNV.all"] <- nrow(x)
        counts["Total_SNV.confident"] <- sum(confident)
       
        # Mutation agreememnt in confident calls
        counts["SNV_agreed_bases.confident"] <- sum(confident & (mutA == mutB | mutA == mutBinverted))
        counts["SNV_agreed_direction.confident"] <- sum(confident & mutA == mutB)
        counts["SNV_disagreed.confident"] <- sum(confident & mutA != mutB & mutA != mutBinverted & !noCallMutB )
        #counts["SNV_present_only_in_first.confident"] <- 

        # Mutation agreement in all calls
        counts["SNV_agreed_bases.all"] <- sum(mutA == mutB | mutA == mutBinverted)
        counts["SNV_agreed_direction.all"] <- sum(mutA == mutB)
        counts["SNV_disagreed.all"] <- sum( mutA != mutB & mutA != mutBinverted & !noCallMutB )
        #counts["SNV_present_only_in_first.all"] <- sum(grepl("[Xx]", mutB))
        
        # Get frequency tables of alt counts and maf for SNV_agreed_direction.all
        if(sum(mutA != mutB) >0) {
            x[ mutA != mutB, ALT_I[2]] <- 0 # convert alt counts of mismatches t "0" for frequency table 
        }
        
        alt_count_freq <- table(x[, ALT_I[2]])
        maf_freq <- table(round(x[, ALT_I[2]] / x[, COV_I[2]] + 5e-3, 2))
    }
    

    return(list(counts, alt_count_freq, maf_freq))
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

simplifyMutationTypes <- function(mutCount) {

    # From a data.frame of counts of SNVs, reduces the number or rows
    # by merging equivalent mutation type counts into one (e.g. C>T and G>A)
    # 
    # @args mutCount - data.frame - table of mutations counts, output of getCountsByType()
    # 
    # @return - data.frame - table of mutations counts with less row number

    results <- list()
    for (i in 1:nrow(MUT_KEY)) {
        mutPair <- MUT_KEY[i,]
        currentMutCount <- mutCount[mutCount$mutType  %in% mutPair, , drop = F]
        currentMutCount <- currentMutCount[, -ncol(currentMutCount)]
        currentMutCount <- colSums(currentMutCount)
        results[[i]] <- currentMutCount
    }

    results <- as.data.frame(do.call(rbind, results))
    results$mutType <- MUT_KEY[,1]

    return(results)
}

main()

# Takes mutation files with two experiments and counts how many times mutations
# are found in both experiments

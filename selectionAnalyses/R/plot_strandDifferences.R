# Plots mutation type counts by strand
# take the mutation table output from dnds script

library("ggplot2")
library("dplyr")
library("tidyr")
source("../../R/ggthemes.R")
source("../../R/gtex.R", chdir = T)

nucleotideFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n1"
NUC_COUNTS <- read.table(nucleotideFile, sep = "\t", stringsAsFactors = F, header = F, row.names = 1)

args <- commandArgs(T)
mutationFile <- args[1] #
outPlot1 <- args[2] # normalized
outPlot2 <- args[3] # raw counts
outPlot3 <- args[4] # normalized - foldchange box
outPlot4 <- args[5] # raw counts - foldchange box
outPlot5 <- args[6] # raw counts - foldchange box
outFile <- args[7] # raw counts - foldchange box
outFile2 <- args[8] # raw counts - foldchange box
outFile3 <- args[9] # strand assymetries per gene
#mutationFile <- c("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Spleen/n6_0.0_0.7/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt")

mutationTypes <- c(`G>A` = "C>T", `G>T` = "C>A", `G>C` = "C>G", `A>G` = "T>C", `A>C` = "T>G", `A>T` = "T>A")


#------------------------
# METHODS

getMutationStats <- function(x, normalized = F) {
    
    mutationStats_all <- x %>%
                         #Counts mutations per person
                         group_by(sampleID, onTranscribed, mutation) %>%
                         count(mutation) 
                         
    if(normalized) {
        mutationStats_all$n <- mutationStats_all$n / NUC_COUNTS[ gsub("(\\w)>.+", "\\1", mutationStats_all$mutation), 1] * 2
    }
    
    mutation_paired_fc <- mutationStats_all %>%
                         #Calculates mean and stderr across all people
                         group_by(sampleID, mutation) %>%
                         spread(onTranscribed, n) %>%
                         mutate(fc = `TRUE`/`FALSE`) %>%
                         ungroup()
                         
    mutationStats <- mutationStats_all %>%
                         #Calculates mean and stderr across all people
                         group_by(onTranscribed, mutation) %>%
                         summarise(avg = mean(n), std = sd(n), stderr = sd(n) / length(n))
                     
    return(list(mutationStats,mutation_paired_fc, mutationStats_all))
                     
                         
}

get_strand_bias_per_gene_pvalue <- function(mutationTable) {
    
    per_gene <- mutationTable %>%
        group_by(gene, mutation, sampleID) %>%
        summarise(transcribed = sum(onTranscribed), non_transcribed = sum(!onTranscribed)) %>%
        ungroup() %>%
        group_by(gene, mutation) %>%
        summarise(wilcox.p = wilcox.test(transcribed, non_transcribed)[["p.value"]], 
                  transcribed = sum(transcribed) + 0.9, non_transcribed = sum(non_transcribed) + 0.9, 
                  total = sum(transcribed, non_transcribed), 
                  n_samples_test = n(),
                  fc = transcribed / non_transcribed) %>%
        ungroup()
    
    per_gene$fdr <- p.adjust(per_gene$wilcox.p, method = "BH")
    return(per_gene)
}

    
plotBars <- function(x) {
    
    ggplot(x, aes(x = mutation, y = avg, fill = onTranscribed)) + 
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = avg - std, ymax = avg + std), position=position_dodge(width=0.9), width = 0.1) +
    theme_grid_y() +
    scale_fill_manual(values = c("#606060","#939393")) + 
    theme(legend.position="top")
    
}

plotBoxes <- function(x) {
    
    maxVal <- log2(x$fc)
    maxVal <- ceiling(max(abs(maxVal[is.finite(maxVal)])))
    
    ggplot(x, aes(x = mutation, y = log2(fc))) + 
        geom_violin(fill = "grey50") +
        geom_boxplot(width = 0.15, colour = "grey60", fill = "grey80") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        ylab("Fold change mutations (Transcribed/Non-transcribed)") +
        #ylim(-maxVal,maxVal) + 
        theme_grid_y() 
    
}

plotBoxes2 <- function(x) {
    
    ggplot(x, aes(x = mutation, y = n, fill = onTranscribed)) + 
        geom_boxplot() +
        scale_fill_manual(values = c("#606060","#939393")) + 
        ylab("Mutations") +
        #ylim(-maxVal,maxVal) + 
        theme_grid_y() 
    
}


get_age_cor <- function(fc_table) {
    
    
    # arriging and getting ids
    fc_table <- mutation_stats_nonNorm[[2]]
    fc_table <- fc_table[!is.na(fc_table$fc),]
    gtexIds <- sraToGtex(unique(fc_table$sampleID), formatOut = "short")
    gtexIds <- gtexIds[!is.na(gtexIds)]
    
    # reading metadata 
    meta <- readMetadataGTEX("AGE")[gtexIds,,drop =F]
    shared_inds <- gtexIds[gtexIds %in% rownames(meta)]
    meta <- meta[shared_inds,,drop = F]
    
    # merging
    fc_table <- fc_table[fc_table$sampleID %in% names(shared_inds),]
    fc_table$gtexId <- shared_inds[fc_table$sampleID]
    fc_table$age <- meta[fc_table$gtexId, "AGE"]
    
    results <- list()
    for(mut in unique(fc_table$mutation)) {
        current_table <- fc_table[fc_table$mutation == mut, ]
        cor_result <- cor.test(current_table$fc, current_table$age, method = "sp")
        results [[mut]] <- data.frame(mutation = mut, cor = cor_result[["estimate"]], pvalue = cor_result[["p.value"]])
    }
    
    results <- do.call(rbind, results)
    
}



#------------------------
# MAIN

# Reads in mutations
mutationTable <- read.table(mutationFile, sep = "\t", header = T, stringsAsFactors = F)
mutationTable$mutation_cod <- paste0(mutationTable$ref_cod, ">", mutationTable$mut_cod)


# Appending if mutation is on transcribed and creating 6-type mutation profile
mutationTable$onTranscribed <- T
mutationTable[mutationTable$mutation_cod %in% mutationTypes[1:6], "onTranscribed"] <- F

mutationTable$mutation <- mutationTable$mutation_cod
mutationTable[mutationTable$onTranscribed, "mutation"] <- mutationTypes[mutationTable[mutationTable$onTranscribed, "mutation"]]

# Calculating values per gene
per_gene_assymetry <- get_strand_bias_per_gene_pvalue(mutationTable)
per_gene_assymetry <- per_gene_assymetry[order(per_gene_assymetry$mutation, -per_gene_assymetry$fc),]


# Calculating total mutations per person
mutation_stats_norm <-  getMutationStats(mutationTable, normalized = T)
mutation_stats_nonNorm <-  getMutationStats(mutationTable, normalized = F)

# Getting correlation with age 
age_cor <- get_age_cor(mutation_stats_nonNorm[[2]])

p1 <- plotBars( mutation_stats_norm[[1]])
p2 <- plotBars( mutation_stats_nonNorm[[1]])

p3 <- plotBoxes( mutation_stats_norm[[2]])
p4 <- plotBoxes( mutation_stats_nonNorm[[2]])

p5 <- plotBoxes2( mutation_stats_nonNorm[[3]])

ggsave(outPlot1, p1)
ggsave(outPlot2, p2)
ggsave(outPlot3, p3)
ggsave(outPlot4, p4)
ggsave(outPlot5, p5, width = 4)
write.table(mutation_stats_nonNorm[[1]], outFile, sep = "\t", quote = F, row.names = F, col.names = T)
write.table(age_cor, outFile2, sep = "\t", quote = F, row.names = F, col.names = T)
write.table(per_gene_assymetry, outFile3, sep = "\t", quote = F, row.names = F, col.names = T)

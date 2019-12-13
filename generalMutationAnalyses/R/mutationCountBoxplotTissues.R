library("ggplot2")
library("dplyr")
library("ggrepel")
library("metap")
source("../../R/ggthemes.R")
source("../../R/misc.R")
source("../../R/plots.R")

args <- commandArgs(T)

workingDir <- args[1] # Dir with counts, e.g mutationCount/count/
maf <- args[2] # Maf of interest 
out_prefix <- args[3]
read_counts_file <- args[4]
stem_cell_division_file <- args[5]

#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/"
#maf <- "n6_0.0_0.7"
#out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/totalCountsBox/raw_n6_0.0_0.7_"
#out_prefix <- "~/raw_n6_0.0_0.7_"
#read_counts_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"
#stem_cell_division_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/stem_cell_divisions.txt"

WIDTH_PER_BOXPLOT <- 0.2 # In inches


# Read read counts per sample
read_counts <- read.table(read_counts_file, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
# Read stem cell division per sample
stem_cell <- read.table(stem_cell_division_file, sep = "\t", header = T, stringsAsFactors = F)
stem_cell <- stem_cell[!is.na(stem_cell$cancer_type),]
stem_cell$log10_stem_cell_division <- log10(stem_cell$n_divisions_all_stem_cells_per_lifetime)
stem_cell$log10_life_time_cancer_risk <- log10(stem_cell$life_time_cancer_risk)

if(!dir.exists(workingDir))
    stop(paste(workingDir, " does not exist"))

tissues <- basename(list.dirs(workingDir, recursive = F))
if(length(tissues) == 0)
    stop(paste(workingDir, " has no tissue folders"))

# Reading data
allTissueCounts <- list()
for(tissue in tissues) {
    
    if(grepl("EXO", tissue))
        next
    
    tissuePath <- file.path(workingDir, tissue, maf)
    if(!dir.exists(tissuePath))
        next
    
    samples = list.files(tissuePath)
    if(length(tissues) == 0)
        stop(paste(tissuePat, " has no sample count files"))
    
    
    countTable <- data.frame(sample = basename(samples), mutations = 0,
                             T_A = 0, T_G = 0, T_C = 0,
                             C_A = 0, C_T = 0, C_G = 0,
                             stringsAsFactors = F)
    
    for (i in 1:length(samples)) {
        samplePath <- file.path(tissuePath, samples[i])
        mutation <- read.table(samplePath, sep = "\t", stringsAsFactors = F, header = F)
        countTable[i, "sample"] <- gsub("\\.txt", "", samples[i])
        countTable[i, "tissue"] <- tissue
        
        countTable[i, "mutations"] <- sum(mutation)
        
        countTable[i, "T_A"] <- mutation[1,1]
        countTable[i, "T_G"] <- mutation[1,3]
        countTable[i, "T_C"] <- mutation[1,4]
        
        countTable[i, "C_A"] <- mutation[2,1]
        countTable[i, "C_T"] <- mutation[2,2]
        countTable[i, "C_G"] <- mutation[2,3]
        
    }
    
    allTissueCounts[[tissue]] <- countTable
} 

# Merging results
allTissueCounts <- (do.call(rbind, allTissueCounts))
allTissueCounts$read_depth <- 0
allTissueCounts[,"read_depth"] <- read_counts[allTissueCounts$sample, "n_uniqueMapped"]
allTissueCounts <- allTissueCounts[!is.na(allTissueCounts$read_depth), ]

# Getting total counts per tissue 
fun <- function(x) sum(as.numeric(x))
stat_tissue <- allTissueCounts %>%
   group_by(tissue) %>%
   summarise_at(vars(mutations:read_depth), fun) %>%
   ungroup()

#Order by mean mutations
pvalues_stem <- list()
correlation_mut_seqDepth <- list()
for(i in c("mutations", "T_A", "T_G", "T_C", "C_A", "C_T", "C_G")) {

    medianMutations <-sort(tapply(allTissueCounts[, i], allTissueCounts$tissue, median))
    allTissueCounts$tissue <- factor(allTissueCounts$tissue, levels = names(medianMutations), ordered = T)

    p <- ggplot(allTissueCounts, aes_string( y = i, x = "tissue")) +
            geom_boxplot() + 
            scale_y_log10() +
            coord_flip() +
            xlab(i) +
            theme_grid_x()
        
    # Save correlation
    flush.console()
    print(stat_tissue[,"read_depth"])
    print(stat_tissue[,i])
    
    current_cor <- cor.test(stat_tissue[,"read_depth", drop = T], stat_tissue[,i, drop = T], method = "sp")
    label_cor <- paste0("r = ", signif(current_cor[["estimate"]], 2), "\np = ", signif(current_cor[["p.value"]], 2))
    correlation_mut_seqDepth[[i]] <- data.frame(mut_type = i, spearman = current_cor[["estimate"]], pvalue = current_cor[["p.value"]])

    # Plot total number of mutations vs total seq depth
    p2 <- ggplot(stat_tissue, aes_string(x = "read_depth", y = i)) + 
             stat_smooth(method = "lm", se = F, color = "lightskyblue3") +
             #geom_abline(linetype = "dashed", colour = "grey40", intercept = 0, slope = max(stat_tissue[,i])/max(stat_tissue$read_depth)) +
             geom_text_repel(aes(label = tissue), colour = "grey30", fontface = "bold") + 
             annotate("text", x = -Inf, y = Inf, label = label_cor, hjust = 0, vjust = 1) +
             geom_point() + 
             ylab(i) +
             xlab("Sequencing depth") +
             theme_bw()
    

    #Plot predicted vs expected mutations
    predicted <- predict.lm(lm(as.formula(paste(i, "~ read_depth")), data = stat_tissue))
    obs_pred <- data.frame(tissue = stat_tissue$tissue, observed = stat_tissue[,i, drop = T], predicted = predicted)
    obs_pred$ratio <- log2(obs_pred$observed / obs_pred$predicted)
    
    stem_cell_merged <- merge(obs_pred, stem_cell, by = "tissue")
    stem_cell_merged$tissue <- as.character(stem_cell_merged$tissue)
    
    # Add brain
    stem_cell_merged <- rbind(stem_cell_merged, stem_cell_merged[stem_cell_merged$tissue == "Brain_Cortex",])
    stem_cell_merged[nrow(stem_cell_merged), "tissue"] <- "Brain average"
    stem_cell_merged[nrow(stem_cell_merged), "ratio"] <- mean(stem_cell_merged[grep("Brain_", stem_cell_merged$tissue), "ratio"])
    stem_cell_merged <- stem_cell_merged[!grepl("Brain_", stem_cell_merged$tissue),]
    
    # Add Esophagus
    stem_cell_merged <- rbind(stem_cell_merged, stem_cell_merged[stem_cell_merged$tissue == "Esophagus_Mucosa",])
    stem_cell_merged[nrow(stem_cell_merged), "tissue"] <- "Esophagus average"
    stem_cell_merged[nrow(stem_cell_merged), "ratio"] <- mean(stem_cell_merged[grep("Esophagus_", stem_cell_merged$tissue), "ratio"])
    stem_cell_merged <- stem_cell_merged[!grepl("Esophagus_", stem_cell_merged$tissue),]
    
    # Remove esophagus muscularis
    stem_cell_merged <- stem_cell_merged[!grepl("muscular", stem_cell_merged$tissue),]
    
    # Add pair label
    stem_cell_merged$pair <- paste0(stem_cell_merged$tissue, "-", stem_cell_merged$cancer_type)
    
    # Remove pancreas
    #stem_cell_merged <- stem_cell_merged[!grepl("ancre", stem_cell_merged$tissue),]
    
    x1 <- min(stem_cell_merged$log10_stem_cell_division)
    x2 <- max(stem_cell_merged$log10_stem_cell_division)
    y1 <- min(stem_cell_merged$ratio)
    y2 <- max(stem_cell_merged$ratio)
    slo <- slope(x1, x2, y1, y2)
                                                               
    p_stem <- scatter(stem_cell_merged, x = "log10_stem_cell_division", y = "ratio", regression = T, method_cor = "sp")
    p_stem <- p_stem + #geom_abline(slope = slo, intercept = y_intercept(x1, y1, slo), linetype = "dashed", colour = "grey40", ) +
        geom_text_repel(aes(label = tissue)) +
        ylab("Mutations log2(observed / expected)") +
        xlab("log10(Number of stem cell divisions)") +
        theme_noGrid()
    
    p_risk <- scatter(stem_cell_merged, x = "log10_life_time_cancer_risk", y = "ratio", regression = T, method_cor = "sp")
    p_risk <- p_risk + #geom_abline(slope = slo, intercept = y_intercept(x1, y1, slo), linetype = "dashed", colour = "grey40", ) +
        geom_text_repel(aes(label = pair)) +
        ylab("Mutations log2(observed / expected)") +
        xlab("Life time cancer risk") +
        theme_noGrid()
    
    
    stem_cor <- cor.test(stem_cell_merged$log10_stem_cell_division, stem_cell_merged$ratio, method = "sp")
    risk_cor <- cor.test(stem_cell_merged$log10_life_time_cancer_risk, stem_cell_merged$ratio, method = "sp")
    pvalues_stem[[i]] <- data.frame(mut_type = i, 
                                    spearmean_cor_stem = stem_cor[["estimate"]], pvalue_stem = stem_cor[["p.value"]], 
                                    spearmean_cor_risk = risk_cor[["estimate"]], pvalue_risk = risk_cor[["p.value"]],
                                    stringsAsFactors = F)
        

    ggsave(paste0(out_prefix, "_", i, "_box.pdf"), p, height = WIDTH_PER_BOXPLOT * length(unique(allTissueCounts$tissue)))
    ggsave(paste0(out_prefix, "_", i, "_scatter_total.pdf"), p2)
    ggsave(paste0(out_prefix, "_", i, "_stem_cell_divisions.pdf"), p_stem)
    ggsave(paste0(out_prefix, "_", i, "_stem_cell_cancer_risk.pdf"), p_risk)
}


# Gathering correlations seq depth and mutations
correlation_mut_seqDepth <- do.call(rbind, correlation_mut_seqDepth)

# Gattering pvalues for stem cell correlations
pvalues_stem <- do.call(rbind, pvalues_stem)
pvalues_stem$p_fisher_stem <- sumlog(pvalues_stem[pvalues_stem$mut_type != "mutations", "pvalue_stem"])$p
pvalues_stem$p_fisher_risk <- sumlog(pvalues_stem[pvalues_stem$mut_type != "mutations", "pvalue_risk"])$p

write.table(pvalues_stem, paste0(out_prefix, "_stem_correlations_pvalues.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(correlation_mut_seqDepth, paste0(out_prefix, "_seqDepth_correlations_pvalues.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

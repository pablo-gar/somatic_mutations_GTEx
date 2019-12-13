library("dplyr")

args <- commandArgs(T)

workingDir <- args[1] # Dir with counts, e.g mutationCount/count/
maf <- args[2] # Maf of interest 
read_counts_file <- args[3]
out_file <- args[4]

#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/"
#maf <- "n6_0.0_0.7"
#out_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/observed_predicted_counts_per_tissue/predicted_observed.txt"
#read_counts_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"

WIDTH_PER_BOXPLOT <- 0.2 # In inches


# Read read counts per sample
read_counts <- read.table(read_counts_file, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

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


# Getting predicted and observed values 
results <- list()
for(i in c("mutations", "T_A", "T_G", "T_C", "C_A", "C_T", "C_G")) {


    #Plot predicted vs expected mutations
    predicted <- predict.lm(lm(as.formula(paste(i, "~ read_depth")), data = stat_tissue))
    obs_pred <- data.frame(tissue = stat_tissue$tissue, mut_type = i, reads = stat_tissue$read_depth, observed = stat_tissue[,i, drop = T], predicted = predicted)
    obs_pred$ratio <- obs_pred$observed / obs_pred$predicted
    
    results[[i]] <- obs_pred
}

results <- do.call(rbind, results)
write.table(results, out_file, sep = "\t", quote = F, row.names = F, col.names = T)

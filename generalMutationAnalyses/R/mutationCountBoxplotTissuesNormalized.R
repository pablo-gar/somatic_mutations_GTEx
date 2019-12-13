library(ggplot2)
source("../../R/ggthemes.R")

args <- commandArgs(T)

workingDir <- args[1] # Dir with counts, e.g mutationCount/count/
maf <- args[2] # Maf of interest 
readCountFile <- args[3]
outplot <- args[4]

WIDTH_PER_BOXPLOT <- 0.5 # In inches

if(!dir.exists(workingDir))
    stop(paste(workingDir, " does not exist"))


tissues <- basename(list.dirs(workingDir, recursive = F))
if(length(tissues) == 0)
    stop(paste(workingDir, " has no tissue folders"))

# Reading data

# Read read counts
readCountTable <- read.table(readCountFile, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

# Read mutation counts 
allTissueCounts <- list()
for(tissue in tissues) {
    
    if(grepl("EXO", tissue))
        next
    
    tissuePath <- file.path(workingDir, tissue, maf)
    if(!dir.exists(tissuePath))
        next
    
    if(grepl("Whole_Blood_EXO", tissuePath))
        next
    
    samples = list.files(tissuePath)
    if(length(tissues) == 0)
        stop(paste(tissuePat, " has no sample count files"))
    
    
    countTable <- data.frame(sample = basename(samples), mutations = 0, n_reads = 0, stringsAsFactors = F)
    for (i in 1:length(samples)) {
        
        
        sampleId <- gsub("\\..*", "",  basename(samples[i]))
        
        if(sampleId %in% rownames(readCountTable)) {
        
            samplePath <- file.path(tissuePath, samples[i])
            mutation <- read.table(samplePath, sep = "\t", stringsAsFactors = F, header = F)
            countTable[i, "sample"] <- samples[i]
            countTable[i, "mutations"] <- sum(mutation)
            countTable[i, "tissue"] <- tissue
            countTable[i, "n_reads"] <- readCountTable[sampleId, "n_reads"]
        }
                
    }
    
    
    # Assesing if we have mutations and read number
    countTable <- countTable[countTable$n_reads != 0,]
    if(nrow(countTable) > 2)
        allTissueCounts[[tissue]] <- countTable
}

allTissueCounts <- (do.call(rbind, allTissueCounts))
allTissueCounts$mutations = allTissueCounts$mutations/ allTissueCounts$n_reads

#Order by mean mutations
medianMutations <-sort(tapply(allTissueCounts$mutations, allTissueCounts$tissue, median))
allTissueCounts$tissue <- factor(allTissueCounts$tissue, levels = names(medianMutations), ordered = T)

p <- ggplot(allTissueCounts, aes( y = mutations/n_reads, x = tissue)) +
geom_boxplot() +
ylab("Mutations per mapped read") +
xlab("Tissue") +
coord_flip() +
theme_grid_x()

ggsave(outplot, p, height = WIDTH_PER_BOXPLOT * length(unique(allTissueCounts$tissue)))

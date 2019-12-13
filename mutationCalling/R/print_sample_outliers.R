# Takes a folder and traverses all subfolders until it finds all SRR*txt files
# the initial folder is supposed to contain all the mutation table counts (2x4 mutation matrix)
# 
# Prints to screen outliers with very low or high number of mutations based on these rules
# Save a table with information of outlier samples
# 
# - minimum number of mutations: 50
# - max number of mutations: TODO
#
# Assumes a that the dir is ./tissue/maf/sample.txt
# Usage Rscript print_sample_outliers.R /scratch/users/paedugar/somaticMutationsProject/mutationCount/count/ outTable.txt

min_mutations <- 15 
max_mutations <- 7500

cmdArgs <- commandArgs(T)

workingDir <- cmdArgs[1] # Dir with counts, e.g mutationCount/count/
outTable <- cmdArgs[2]

#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count"
#outTable <- "/scratch/users/paedugar/somaticMutationsProject/purgedSamples.txt"

if(!dir.exists(workingDir))
    stop(paste(workingDir, " does not exist"))

samples <- file.path(workingDir, list.files(workingDir, recursive = T))
samples <- samples[!grepl("EXO",samples)]
    
countTable <- data.frame(sample = basename(samples), mutations = 0, tissue = "", maf = "", stringsAsFactors = F)
for (i in 1:length(samples)) {
    mutation <- read.table(samples[i], sep = "\t", stringsAsFactors = F, header = F)
    countTable[i, "sample"] <- basename(samples[i])
    countTable[i, "mutations"] <- sum(mutation)
    countTable[i, "tissue"] <- basename(dirname(dirname(samples)))[i]
    countTable[i, "maf"] <- basename(dirname(samples))[i]
}



# Take top and bottom samples and print them
bottom <- countTable[ countTable$mutations < min_mutations, , drop = F]
top <- countTable[ countTable$mutations > max_mutations, , drop = F]

if(nrow(bottom) > 0)
    bottom$type <- "bottom"
if(nrow(top) > 0)
    top$type <- "top"

allPurge <- rbind(top,bottom)

write.table(allPurge, outTable, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
cat(paste(allPurge$sample, collapse = " "))

# Removes files that did not pass 
mutation_root_dir <- dirname(workingDir)
allFiles <- file.path(mutation_root_dir, list.files(mutation_root_dir, recursive = T))

for(i in seq_len(nrow(allPurge))) {
    
    prefix <- file.path(allPurge[i, "maf"], allPurge[i, "sample"])
    to_remove <- allFiles[grep(prefix, allFiles)]
    cat(paste(to_remove, collapse = "\n"), "\n")
    
    file.remove(to_remove)
    
}

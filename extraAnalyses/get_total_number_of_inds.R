# Prints number of individuals

source("../R/gtex.R", chdir = T)

# Reading files
mutation_count_dir <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count"
mutation_files <- list.files(mutation_count_dir, full = T, recursive = T)
mutation_files <- mutation_files[!grepl("EXO", mutation_files)]

# Samples
tissues <- unique(basename(dirname(dirname(mutation_files))))
samples <- gsub(".txt", "", basename(mutation_files))
inds <-  sraToGtex(samples, formatOut="short")

cat("N samples = ", length(samples), "\n")
cat("N tissues = ", length(tissues), "\n")
cat("N inds = ", length(unique(inds)), "\n")

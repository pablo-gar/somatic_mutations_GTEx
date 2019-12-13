# Creates a table with the number of samples per tissue
source("../../R/gtex.R", chdir = T)

out_file <- file.path(CONFIG$scratchDir, CONFIG$generalMutationAnalyses$root, "table_tissue_sample_count", "tissue_sample_count.txt")

samples <- list.files(file.path(CONFIG$scratchDir, CONFIG$mutationCountDir$count), recursive = T, full = T, pattern = "*txt")
tissues <- basename(dirname(dirname( samples)))
tissues <- as.data.frame(table(tissues))

dir.create(dirname(out_file), recursive = T)

write.table(tissues, out_file, sep = "\t", quote = F, col.names = F, row.names = F)

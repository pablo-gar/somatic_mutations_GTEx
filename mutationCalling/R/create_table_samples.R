source("../../R/gtex.R", chdir = T)

samples <- list.files(file.path(CONFIG$scratchDir, CONFIG$mutationCountDir$count), recursive = T, full = T)
samples <- gsub(".txt", "", basename(samples))
out_table_file <- file.path(CONFIG$auxiliaryFiles$sampleTable)
readCountFile <- CONFIG$auxiliaryFiles$readCountTable
transDiversityFile <- CONFIG$auxiliaryFiles$transcriptomeDiversity


# Read is info
out_table <- data.frame(sraId = samples,
                        gtexId = sraToGtex(samples, formatOut = "long"),
                        gtexInd = sraToGtex(samples, formatOut = "short"),
                        tissue =  querySraRunTable(samples, refCol = 16 , queryCol = 21),
                        stringsAsFactors = F
                        )

# Read metadata
meta <- readMetadataGTEX(c("AGE", "RACE", "BMI", "GENDER"))
meta <- meta [rownames(meta) %in% out_table$gtexInd,]

out_table <- cbind(out_table, meta[out_table$gtexInd,])

# Read sample anno
meta <- readSampleAnnotationGTEX(c("SMRIN", "SMTSISCH", "SMTSPAX"))

out_table <- cbind(out_table, meta[out_table$gtexId,])

# Read genotypes
genotype <- as.data.frame(t(readGenotypePCA())[,1:3])

out_table <- cbind(out_table, genotype[out_table$gtexInd,])

# Read seq depth
seq_depth <- read.table(readCountFile, sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
seq_depth <- seq_depth[,"n_uniqueMapped", drop = F]

out_table <- cbind(out_table, seq_depth[out_table$sraId,])
colnames(out_table)[ncol(out_table)] <- "n_uniqueMapped"

# Read trans diversity 
trans_diver <- read.table(transDiversityFile, sep = "\t", stringsAsFactors = F, header = T, row.names = 1)

out_table <- cbind(out_table, trans_diver[out_table$gtexId,])
colnames(out_table)[ncol(out_table)] <- "transcriptome_diversity"

# Fixing tissue name
out_table$tissue <- gsub("-", "", out_table$tissue)
out_table$tissue <- gsub("\\(", "", out_table$tissue)
out_table$tissue <- gsub("\\)", "", out_table$tissue)
out_table$tissue <- gsub("\\s+", "_", out_table$tissue)

write.table(out_table, out_table_file, sep = "\t", quote = F, col.names = T, row.names = F)

# Makes table of blood transcriptomes for https://cibersort.stanford.edu/runcibersort.php
# Rscript make_table_CIBERSORT.R "Whole Blood" /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/TPM_whole_blood.txt

source("../R/gtex.R")
source("../R/geneTools.R")

cmdArgs <- commandArgs(T)
tissue_interest <- cmdArgs[1]
out_table <- cmdArgs[2]

#tissue_interest <- "Whole Blood"
#out_table <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/TPM_Whole_Blood_allGenes.txt"

#tissue_interest <- "Lung"
#out_table <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/TPM_Lung_allGenes.txt"

#tissue_interest <- "Spleen"
#out_table <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/TPM_Spleen_allGenes.txt"

sraSamples_file <- "/scratch/PI/hbfraser/gtex/anno/GtexSraRunTable.txt"
signatureFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/LM22.txt"

# Get blood samples
samples <- read.table(sraSamples_file, sep = "\t", stringsAsFactors = F, header = T)
samples <- samples[,c("Run_s", "biospecimen_repository_sample_id_s", "body_site_s", "Assay_Type_s")]
samples <- samples[samples$Assay_Type_s == "RNA-Seq" & samples$body_site_s == tissue_interest,]

# Get rna-seq samples
exp_mat <- readAllGtexExpression(samples$biospecimen_repository_sample_id_s)
exp_mat[,1] <- gsub("\\.\\d+","",exp_mat[,1])
hugo_symbols <- ensemblToAlias(exp_mat[,1])
exp_mat[,1] <- hugo_symbols
exp_mat <- exp_mat[!is.na(exp_mat)[,1],]

# Get signatures files
signatures <- read.table(signatureFile, sep = "\t", stringsAsFactors = F, header = T)


# Get overlaping genes
colnames(exp_mat)[1] <- colnames(signatures)[1]
#exp_mat <- exp_mat[exp_mat[,1] %in% signatures[,1], ]
exp_mat <- exp_mat[exp_mat[,1] != "", ]
exp_mat <- exp_mat[!is.na(exp_mat[,1]), ]
signatures <- signatures[signatures[,1] %in% exp_mat[,1], ]

write.table(exp_mat, out_table, sep = "\t", quote = F, row.names = F, col.names = T)
write.table(signatures, paste0(signatureFile, ".", tissue_interest, ".ready"), sep = "\t", quote = F, row.names = F, col.names = T)



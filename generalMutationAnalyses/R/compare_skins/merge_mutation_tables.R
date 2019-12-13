source("../../../R/gtex.R", chdir = T)
library("dplyr")
library("tidyr")

cmdArgs <- commandArgs(T)

mut_table_file_1 <- cmdArgs[1]
mut_table_file_2 <- cmdArgs[2]
out_table_file <- cmdArgs[3]

mut_table_file_1 <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Skin_Sun_Exposed_Lower_leg-n6_0.0_0.7.txt"
mut_table_file_2 <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Skin_Not_Sun_Exposed_Suprapubic-n6_0.0_0.7.txt"
out_table_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates_skin_fc/skin_C_T_fc_covariates.txt"



# Read tables
mut_table_1 <- read.table(mut_table_file_1, sep = "\t", stringsAsFactors = F, header = T)
mut_table_2 <- read.table(mut_table_file_2, sep = "\t", stringsAsFactors = F, header = T)

mut_table_1$tissue <- gsub("-.+", "", basename(mut_table_file_1))
mut_table_2$tissue <- gsub("-.+", "", basename(mut_table_file_2))

# Join and calculate C>T freq and fc
mut_table <- rbind(mut_table_1, mut_table_2)
gtexIds <- sraToGtex(mut_table$sample)
mut_table$gtexId <- gtexIds[mut_table$sample]
mut_table$C_T.freq <- mut_table$C_T / mut_table$all

mut_table <- mut_table[,c("AGE", "RACE", "GENDER", "BMI", "C1", "C2", "C3", "tissue", "gtexId", "C_T.freq")]
mut_table <- mut_table %>%
    spread(tissue, C_T.freq) %>%
    filter(!is.na(Skin_Not_Sun_Exposed_Suprapubic) & !is.na(Skin_Sun_Exposed_Lower_leg)) %>%
    mutate(fc = Skin_Sun_Exposed_Lower_leg / Skin_Not_Sun_Exposed_Suprapubic)

write.table(mut_table, out_table_file, sep = "\t", quote = F, col.names = T, row.names = F)

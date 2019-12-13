
source("../../../R/mutationGeneAnnotation.R", chdir = T)
source("../../../R/geneTools.R", chdir = T)
source("../../../R/gtex.R", chdir = T)

cmdArgs <- commandArgs(T)
output_table <- cmdArgs[1]
output_table_expression <- cmdArgs[2]
driverMutationsFile <- cmdArgs[3]
mutationFiles <- cmdArgs[-(1:3)]


#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.7/"
#mutationFiles <- file.path(workingDir, list.files(workingDir))
#driverMutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/tcga_driverMutations.txt"


#----------------------------------#
# MAIN
#----------------------------------#

####
# Read driver mutations
driverMutations <- read.delim(driverMutationsFile, sep = "\t", stringsAsFactors = F, header = T)
driverMutations$gene_id <- queryRefCDS(driverMutations$Gene, referenceField = "gene_name", queryField = "gene_id")

#if(all(colnames(driverMutations) != "Mutation"))
#    driverMutations <- "no_driver_mutations"

####
# Read mutations and process mutations
mutationsAll <- list()
for(i in 1:length(mutationFiles)) {
    
    mutationFile <- mutationFiles[i]
    
    currentSample <- gsub(".txt", "", basename(mutationFile))
    mutations <- read.table(mutationFile, sep = "\t", stringsAsFactors = F)
    mutations$sample <- currentSample
    
    mutationsAll[[i]] <- mutations
    
}

mutationsAll <- do.call(rbind, mutationsAll)
mutationsAll <- annotateMutations(mutationsAll)

# Select only driver genes
mutationsAll <- mutationsAll[mutationsAll$gene_id %in% driverMutations$gene_id,, drop = F]
driverMutationsAll <- driverMutations
driverMutations <- driverMutations[driverMutations$gene_id %in% mutationsAll$gene_id,, drop = F]

###
# Annotate mutations as driver or not
mutationsAll$driver <- FALSE
for(gene in unique(driverMutations$gene_id)) { 
    current_aa_changes <- mutationsAll[mutationsAll$gene_id %in% gene, "aachange"]
    mutationsAll[mutationsAll$gene_id %in% gene, "driver"] <- current_aa_changes %in% driverMutations[driverMutations$gene_id %in% gene, "Mutation"]
}

####
# Read gene expression
gtex_ids <- sraToGtex (unique(mutationsAll$sample), formatOut = "long")
expression_mat <- readAllGtexExpression(samples = gtex_ids, CONFIG$auxiliaryFiles$readsGenesAllTissues)
rownames(expression_mat) <- gsub("\\..+", "", expression_mat[,1])
expression_mat <- expression_mat[,-1]
expression_mat <- expression_mat[rownames(expression_mat) %in% driverMutationsAll$gene_id,]
expression_stats <- data.frame(gene_id = rownames(expression_mat),
                               gene_name = queryRefCDS(rownames(expression_mat), referenceField = "gene_id", queryField= "gene_name"), 
                               n_ind = ncol(expression_mat),
                               total_reads = apply(expression_mat, 1, sum),
                               median_reads = apply(expression_mat, 1, median),
                               mean_reads = rowMeans(expression_mat),
                               var_reads = apply(expression_mat, 1, var)
                               )


write.table(mutationsAll, output_table, sep = "\t", col.names = T, row.names = F, quote = F)
write.table(expression_stats, output_table_expression, sep = "\t", col.names = T, row.names = F, quote = F)

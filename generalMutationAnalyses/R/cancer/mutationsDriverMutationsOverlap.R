
source("../../../R/mutationGeneAnnotation.R", chdir = T)
source("../../../R/geneTools.R", chdir = T)

cmdArgs <- commandArgs(T)
output_table <- cmdArgs[1]
driverMutationsFile <- cmdArgs[2]
mutationFiles <- cmdArgs[-(1:2)]


#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.5/"
#mutationFiles <- file.path(workingDir, list.files(workingDir))
#driverMutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/tcga_driverMutations.txt"


#----------------------------------#
# MAIN
#----------------------------------#

####
# Read driver mutations
driverMutations <- read.delim(driverMutationsFile, sep = "\t", stringsAsFactors = F, header = T)
driverMutations$gene_id <- aliasToEnsembl(driverMutations$Gene)

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
mutationsAll <- mutationsAll[mutationsAll$gene_id %in% driverMutations$gene_id,]
driverMutations <- driverMutations[driverMutations$gene_id %in% mutationsAll$gene_id,]

###
# Annotate mutations as driver or not
mutationsAll$driver <- FALSE
for(gene in unique(driverMutations$gene_id)) { 
    current_aa_changes <- mutationsAll[mutationsAll$gene_id %in% gene, "aachange"]
    mutationsAll[mutationsAll$gene_id %in% gene, "driver"] <- current_aa_changes %in% driverMutations[driverMutations$gene_id %in% gene, "Mutation"]
}

write.table(mutationsAll)

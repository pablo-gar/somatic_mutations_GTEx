# Performs linear regressions between the expression of all genes and
# and the mutation signature loading for a given cancer type

#library("MatrixEQTL")
library("dplyr")
library("tidyr")
#library("ggplot2")
#library("metap")
#source("../../R/FDR.R")
#source("../../R/plots.R")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)
#source("../../R/mutationGeneAnnotation.R", chdir = T)

cmdArgs <- commandArgs(T)

mutationsFile <- cmdArgs[1]
pathway_gene_file <- cmdArgs[2]
targetGenes <- parseArg(cmdArgs[3], sep = ",", trim = T)
nPerm <- as.numeric(cmdArgs[4])
output <- cmdArgs[5]
output_group_enrichment <- cmdArgs[6]
output_group_perm_tables <- cmdArgs[7]

ignoreCovariates <- c("ETHNCTY")


#pvalueCutoff <- 0.05
#pathway_gene_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/repair_genes.txt"
#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Whole_Blood-n6_0.0_0.7.txt"
#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Spleen-n6_0.0_0.7.txt"
#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Colon_Transverse-n6_0.0_0.7.txt"
#targetGenes <- parseArg("ENSG00000128383,ENSG00000062822,ENSG00000174371,ENSG00000095002,ENSG00000113318,ENSG00000116062,ENSG00000076242,ENSG00000122512,ENSG00000139618,ENSG00000012048,ENSG00000179750,ENSG00000170734,ENSG00000154767,ENSG00000136936,ENSG00000167986,ENSG00000134574,ENSG00000140398,ENSG00000070501,ENSG00000166169,ENSG00000114026,ENSG00000065057,ENSG00000154328,ENSG00000132781", sep = ",", trim = T)
#geneNames <- "APOBEC3A,MLH1,PMS1,BRCA2,BRCA1,APOBEC3B,POLH,XPC,DDB2,NEIL1"
#output <- "/scratch/users/paedugar/somaticMutationsProject/ngeneralMutationAnalyses/results/mutations_expressionAssociation_targeted/Whole_Blood/6_0.0_0.7_all.txt"

#-------------------
# METHODS

rankitNormalize <- function(x, IND = 1) {

    # Normalizes rows (IND = 1) or columns (IND = 2)
    # to a quantile standard normalization (i.e. rankit)

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))

    x <- apply(x, IND, function(x) qnorm((rank(x) - 0.5) / length(x)))
    if(IND == 1)
        x <- t(x)

    return(x)

}

read_H_mat <- function(x) {
    x <- read.table(x, stringsAsFactors = F, header = T, sep = "\t")
    rownames(x) <- 1:nrow(x)
    return(as.matrix(x))
}



#-------------------



#-------------------
# MAIN

# Read mutation mat
mutationMat_withCov <- read.table(mutationsFile, header = T, sep = "\t", stringsAsFactors = F)
mutationMat_withCov <- mutationMat_withCov[,!colnames(mutationMat_withCov) %in% c("C_C", "T_T")]

# Read pathway genes
pathway_gene_table <- read.table(pathway_gene_file, sep = "\t", stringsAsFactors = F,  header = T)
pathway_gene_table$gene <- queryRefCDS(pathway_gene_table$alias, reference = "gene_name", query = "gene_id")

# Eliminate inds for which we don't have all the metadata
mutationMat_withCov <- mutationMat_withCov[rowSums(is.na(mutationMat_withCov)) == 0,]
sample_col <- which(colnames(mutationMat_withCov) ==  "sample")

# Getting ids
sraIds <- mutationMat_withCov$sample
gtexIds <- sraToGtex(sraIds, formatOut = "short")
gtexIdsLong <- sraToGtex(sraIds, formatOut = "long")
rownames(mutationMat_withCov) <- gtexIdsLong

# Select mutation types only
mutationMat <- t(mutationMat_withCov[,1:(sample_col-1)])

#Read metadata
covariates <- t(mutationMat_withCov[(sample_col + 1) : ncol(mutationMat_withCov)])
covariates <- covariates[ !rownames(covariates) %in% ignoreCovariates, , drop = F] 

# Reads expression file
expressionMat <- readAllGtexExpression(gtexIdsLong)
rownames(expressionMat) <- gsub("\\..+", "", expressionMat[,1])

# Getting ids for which we have info for all
sharedInds <- gtexIdsLong %in% colnames(expressionMat) & 
    gtexIdsLong %in% colnames(covariates) 

sraIds <- sraIds[sharedInds]
gtexIdsLong <- gtexIdsLong[sraIds]


# Only genes expressed in >20% of samples
expressionMat <- expressionMat[,gtexIdsLong]
genes_expressed <- rowSums(expressionMat > 0) >= (ncol(expressionMat) * 0.05)
genes_interest <- rownames(expressionMat) %in% targetGenes
expressionMat <- expressionMat[genes_interest | genes_expressed, ]

# SELECT GENES of interest
targetGenes_ind <- rownames(expressionMat) %in% targetGenes
targetGenes_n <- sum(targetGenes_ind)

mutationMat <- mutationMat[,gtexIdsLong]
covariates <- covariates[,gtexIdsLong]

# Eliminate covariates missing or not variable
covariates <- covariates[rowSums(!is.na(covariates)) == ncol(covariates),]
mostRepeated <- apply(covariates,1,function(x) max(table(x)))
covariates <- covariates[mostRepeated/ncol(covariates) < 0.95,]



#' Performs a GWAS between all SNP with MAF >0.1 and clonality
#' mutation signature
#'
#' 
#' The pvalue for the strongest hit is calculated as the frequency of pvalues lower than it in the
#' permutations
#'
#' @param character - path to file containing all the mutations
#' @param character - path to metadata file 
#' @param integer - number of permutations to perform
#' @param character - tissue name
#' @param character - the column name from the mutation file that contains the values to do GWAS on
#' @param character - path to vcf file with a field "GT" with genotypes of form "0/1"
#' @param character - path to output prefix
#'
#' @return two files, the first with suffix `pvals` contains is a table of the strongest hits; the second one with prefix `ntests` contain the number of snp tested per sginature
#'
#' @usage Rscript gwas_clonality.R 


if(!dir.exists("../../R"))
    stop("Can't find the shared R path, make sure to run script from directory where it lives")

library("rSubmitter")
library("dplyr")
library("tidyr")

source("association_functions.R")
source("../../R/gtex.R", chdir = T)
source("../../R/GWAS_methods.R", chdir = T)
source("../../R/misc.R", chdir = T)

MUT_TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "all")

# Gathering arguments
args <- commandArgs(T)
mutation_map_file <- args[1]
metadata_file <- args[2]
vcfFile <-  args[3] 
permutations <- as.numeric(args[4]) 
tissue <- args[5]
phenotype_var <- args[6]
outFile <-  args[7] 


yieldSize <- 2e6
genotypeMode <- "GT"
minorAlleleFreq <- 0.1
topNsignificant <- 100

# For rSubmitter
GWASsource <- "~/scripts/FraserLab/somaticMutationsProject/R/GWAS.R"
extraBashLines <- c("module load mariadb/10.2.11", "module load R/3.4.0")
packages <- "VariantAnnotation"

#mutation_map_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/all_mutations_annotated.txt"
#metadata_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/metadata.txt"
#tissue <- "Prostate"
#permutations <- 1000
#phenotype_var <- "vaf_corrected"
#library("rjson")
#CONFIG <- fromJSON(file = "../../config.json")
#vcfFile <-  CONFIG$vcfFileCommon


#---------------------------
# MAIN

# Reading mutation data and metadata
mutation_map <- read.table(mutation_map_file, header=T, sep="\t", stringsAsFactors=F)

covariates <- read.table(metadata_file, header=T, sep="\t", stringsAsFactors=F)
covariates <- covariates[covariates$RACE == 3,]
rownames(covariates)<-covariates$sraId
covariates <- select(covariates, -sraIds) 

# Re-arraging data
mutation_map <- get_vaf_per_sample(mutation_map, phenotype_var) 
mutation_map$gtexIds <- gtexLongToShort(mutation_map$gtexIds_samples)


# Creating bins of genome
genomeRanges <- binGenome(CONFIG$auxiliaryFiles$hg19_chromSizes, yieldSize, exclude = c("Y", "X"))

################
# Goes over all mutation signatures

allSignatureGWAS <- list()

current_mat <- mutation_map[mutation_map$tissue == tissue,]
current_mat <- current_mat[current_mat$sraIds %in% rownames(covariates),]

current_cov <- t(covariates[current_mat$sraIds,])
current_cov <- current_cov[apply(current_cov, 1, function(x) length(unique(x)) > 1),]
colnames(current_cov) <- current_mat$gtexIds
current_cov <- rankitNormalize(current_cov, IND=1)


allSignatureGWAS <- list()
for (i in MUT_TYPES[MUT_TYPES %in% colnames(mutation_map)]) {
    
    
    phenoVector <- setNames(current_mat[,i,drop=T], current_mat$gtexIds)
    phenoVector <- rankitNormalize_vector(phenoVector)

    # Makes permuted table 
    permuted <- do.call(rbind, lapply(rep(T, permutations), function(x,y) sample(y), y = phenoVector))
    # Real data is first row
    permuted <- rbind(phenoVector, permuted)

    # Doing GWAS
    #permutedResults <- lapply(genomeRanges[1:3], lapplyHelper, vcfFile = vcfFile, phenotype = permuted, genotype = genotypeMode, topNsignificant = topNsignificant, covariateMatrix = covariates, minorAlleleFreq = minorAlleleFreq)
    permutedResults <- rSubmitter::superApply(genomeRanges, lapplyHelper, vcfFile = vcfFile, phenotype = permuted, 
                                  genotype = genotypeMode, topNsignificant = topNsignificant, 
                                  covariateMatrix = current_cov, minorAlleleFreq = minorAlleleFreq,
                                  tasks = 100, mem = "4G", partition = "hbfraser", time = "30:00",
                                  sources = GWASsource, packages = packages, extraBashLines = extraBashLines, 
                                  workingDir = "/scratch/users/paedugar/somaticMutationsProject/tempDir/" )

    mergedResults <- mergeGWASresults(permutedResults, topNsignificant)
    
    mergedResults$realTop$pvalue_perm <- c(sum(mergedResults$permutedMin[-1] <= mergedResults$realTop$pvalue[1] ) / length(mergedResults$permutedMin[-1]), rep(NA, topNsignificant-1))
    mergedResults$realTop$pvalue_bonf <- mergedResults$realTop$pvalue * mergedResults$totalSNPs
    mergedResults$realTop$signature <- i
    mergedResults$realTop$tissue <- tissue
    mergedResults$realTop$nTests <- mergedResults$totalSNPs
    
    allSignatureGWAS[[i]] <- mergedResults
}



pvals <- do.call(rbind, lapply(allSignatureGWAS, function(x) x$realTop))

# Writing results
write.table(pvals, outFile, sep = "\t", row.names = F, quote = F)
#save(allSignatureGWAS, file = "~/allSignatureGWAS.RData")

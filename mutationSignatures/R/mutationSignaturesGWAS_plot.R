# Do module load mariadb/10.2.11
#
#' Performs a GWAS between all SNP with MAF >0.1 and the frequency of the loading for a given 
#' mutation signature
#'
#' Takes the H matrix of loadings -from NMF- and performs the GWAS in each signature individually
#' it performs permutations and captures the lowers pvalue in each permutation
#' 
#' The pvalue for the strongest hit is calculated as the frequency of pvalues lower than it in the
#' permutations
#'
#' @param character - path to H matrix file
#' @param integer - number of permutations to perform
#' @param character - path to vcf file with a field "GT" with genotypes of form "0/1"
#' @param character - path to output prefix
#'
#' @return two files, the first with suffix `pvals` contains is a table of the strongest hits; the second one with prefix `ntests` contain the number of snp tested per sginature
#'
#' @usage Rscript mutationSignaturesGWAS.R 1000 /scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Colon_Sigmoid-n6_0.0_0.5_H.txt /scratch/PI/hbfraser/gtex/raw/Genotype/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_MAF0.01.vcf.gz ~/mutationSignatureGWAS_colon


if(!dir.exists("../../R"))
    stop("Can't find the shared R path, make sure to run script from directory where it lives")

library("VariantAnnotation")
library("ggplot2")

source("../../R/ggthemes.R", chdir = T)
source("../../R/gtex.R", chdir = T)
source("../../R/GWAS_2.R", chdir = T)

# Gathering arguments
args <- commandArgs(T)
Hfile <- args[1] #"/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Prostate-n6_0.0_0.5_H.txt"
permutations <- as.numeric(args[2]) #1000
vcfFile <-  args[3] #CONFIG$vcfFileCommon
outPrefix <-  args[4] #"~/mutationSignatureGWAS"

indInfo <- c("AGE", "BMI", "GENDER")
yieldSize <- 2e6
genotypeMode <- "GT"
minorAlleleFreq <- 0.1
topNsignificant <- 50

# For rSubmitter
GWASsource <- "~/scripts/FraserLab/somaticMutationsProject/R/GWAS.R"
extraBashLines <- c("module load mariadb/10.2.11", "module load R/3.4.0")
packages <- "VariantAnnotation"

#permutations <- 1000
library("rjson")
CONFIG <- fromJSON(file = "../../config.json")
vcfFile <-  CONFIG$vcfFileCommon


Hfile <-"/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Brain_Cerebellum-n6_0.0_0.5_H.txt"
snpId <- "1_76434065_A_G_b37" # cerebellum
signature <- 3
Hfile <-"/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Brain_Frontal_Cortex_BA9-n6_0.0_0.5_H.txt"
snpId <- "8_125433536_G_A_b37"
signature <- 2
Hfile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Prostate-n6_0.0_0.5_H.txt"
snpId <- "9_2803726_C_T_b37" # cerebellum
signature <- 3
Hfile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Adipose_Subcutaneous-n6_0.0_0.5_H.txt"
snpId <- "7_158531177_C_CT_b37" # adipose
signature <- 1

#---------------------------
# METHODS


#---------------------------


#---------------------------
# MAIN

#####
# Reading mutation data
# Reading table
Hmat <- read.table(Hfile, header = T, sep = "\t", stringsAsFactors = T)

# Calculate frequencies
Hmat <- t(t(Hmat) / (colSums(Hmat)))

# Getting ids
sraIds <- colnames(Hmat)
gtexIds <- sraToGtex(sraIds, formatOut = "short")
colnames(Hmat) <- gtexIds
#
##Read metadata
#indMeta <- readMetadataGTEX(indInfo)
#genotypePCA <- t(readGenotypePCA())
#
## Getting ids for which we have info for all
#allInfo <- gtexIds %in% rownames(indMeta) & gtexIds %in% rownames(genotypePCA)
#
#gtexIds<- gtexIds[allInfo]
#
#Hmat <- Hmat[,gtexIds]
#indMeta <- indMeta[gtexIds,]
#genotypePCA <- genotypePCA[gtexIds,1:10]
#covariates <- t(cbind(indMeta, genotypePCA))



phenoVector <- Hmat[signature,]

## Gets individuals in the vcf file 
individualsLogic <- checkIndividuals(vcfFile, findNames(vcfFile, names(phenoVector)))
phenoVector <- phenoVector[individualsLogic]

# Gets genotypes of SNPs
snpGeno <- getSNPsByID (ids = snpId, vcfFile = vcfFile, individuals = names(phenoVector), genotype = genotypeMode, chromFormat = "ENSEMBL")
for(i in 1:ncol(snpGeno)) snpGeno[,i] <- as.character(snpGeno[,i])
colnames(snpGeno) <- paste0("snp_", colnames(snpGeno))
rownames(snpGeno) <- GTEX$simplifyId(rownames(snpGeno))
snpGeno$signature_weight <-  phenoVector[rownames(snpGeno)]


snpGeno <- snpGeno[rowSums(is.na(snpGeno)) ==  0 ,]
ggplot(snpGeno, aes_string(x = paste0("snp_", snpId), y = "signature_weight")) +
geom_boxplot() +
theme_grid_y()





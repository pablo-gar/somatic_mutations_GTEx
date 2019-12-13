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

library("rSubmitter")

source("../../R/gtex.R", chdir = T)
source("../../R/GWAS_methods.R", chdir = T)

# Gathering arguments
args <- commandArgs(T)
Hfile <- args[1] #"/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Liver-n6_0.0_0.5_H.txt"
permutations <- as.numeric(args[2]) #1000
vcfFile <-  args[3] #CONFIG$vcfFileCommon
outPrefix <-  args[4] #"~/mutationSignatureGWAS"

fromConsensus <- T

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
#Hfile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Liver-n6_0.0_0.7_H_consensus.txt"
#library("rjson")
#CONFIG <- fromJSON(file = "../../config.json")
#vcfFile <-  CONFIG$vcfFileCommon


#---------------------------
# MAIN

#####
# Reading mutation data
# Reading table
Hmat <- read.table(Hfile, header = T, sep = "\t", stringsAsFactors = T)
if(fromConsensus) {
    rownames(Hmat) <- Hmat[,1]
    Hmat <- Hmat[,-1]
}

# Calculate frequencies
Hmat <- t(t(Hmat) / (colSums(Hmat)))

# Getting ids
sraIds <- colnames(Hmat)
gtexIds <- sraToGtex(sraIds, formatOut = "short")
colnames(Hmat) <- gtexIds

#Read metadata
indMeta <- readMetadataGTEX(indInfo)
genotypePCA <- t(readGenotypePCA())

# Getting ids for which we have info for all
allInfo <- gtexIds %in% rownames(indMeta) & gtexIds %in% rownames(genotypePCA)

gtexIds<- gtexIds[allInfo]

Hmat <- Hmat[,gtexIds]
indMeta <- indMeta[gtexIds,]
genotypePCA <- genotypePCA[gtexIds,1:3]
covariates <- t(cbind(indMeta, genotypePCA))

        

## Creating bins of genome
genomeRanges <- binGenome(CONFIG$auxiliaryFiles$hg19_chromSizes, yieldSize, exclude = c("Y", "X"))

################
# Goes over all mutation signatures

allSignatureGWAS <- list()

for (i in 1:nrow(Hmat)) {
    phenoVector <- Hmat[i,]

    # Makes permuted table 
    permuted <- do.call(rbind, lapply(rep(T, permutations), function(x,y) sample(y), y = phenoVector))
    # Real data is first row
    permuted <- rbind(phenoVector, permuted)

    # Doing GWAS
    #permutedResults <- lapply(genomeRanges[1:3], lapplyHelper, vcfFile = vcfFile, phenotype = permuted, genotype = genotypeMode, topNsignificant = topNsignificant, covariateMatrix = covariates, minorAlleleFreq = minorAlleleFreq)
    permutedResults <- rSubmitter::superApply(genomeRanges, lapplyHelper, vcfFile = vcfFile, phenotype = permuted, 
                                  genotype = genotypeMode, topNsignificant = topNsignificant, 
                                  covariateMatrix = covariates, minorAlleleFreq = minorAlleleFreq,
                                  tasks = 200, mem = "4G", partition = "hbfraser", time = "30:00",
                                  sources = GWASsource, packages = packages, extraBashLines = extraBashLines, 
                                  workingDir = "/scratch/users/paedugar/somaticMutationsProject/tempDir/" )

    mergedResults <- mergeGWASresults(permutedResults, topNsignificant)
    
    mergedResults$realTop$pvalue_perm <- c(sum(mergedResults$permutedMin[-1] <= mergedResults$realTop$pvalue[1] ) / length(mergedResults$permutedMin[-1]), rep(NA, topNsignificant-1))
    mergedResults$realTop$pvalue_bonf <- mergedResults$realTop$pvalue * mergedResults$totalSNPs
    mergedResults$realTop$signature <- i
    mergedResults$realTop$nTests <- mergedResults$totalSNPs
    
    allSignatureGWAS[[i]] <- mergedResults
}


pvals <- do.call(rbind, lapply(allSignatureGWAS, function(x) x$realTop))

# Writing results
write.table(pvals, paste0(outPrefix, "_pvals.txt"), sep = "\t", row.names = F, quote = F)
#save(allSignatureGWAS, file = "~/allSignatureGWAS.RData")

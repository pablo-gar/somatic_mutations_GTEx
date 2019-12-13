# ONLY USE WITH COVARIATES
# 
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
#' @usage Rscript mutations_GWAS.R /scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates_skin_fc/skin_C_T_fc_covariates.txt 1000 /scratch/PI/hbfraser/gtex/raw/Genotype/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_MAF0.01.vcf.gz /scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_skins/Skins-n6_0.0_0.7.txt


if(!dir.exists("../../R"))
    stop("Can't find the shared R path, make sure to run script from directory where it lives")

library("rSubmitter")
source("../../../R/gtex.R", chdir = T)
source("../../../R/GWAS_methods.R", chdir = T)

# Gathering arguments
args <- commandArgs(T)
mutationsFile <- args[1] 
permutations <- as.numeric(args[2]) #1000
vcfFile <-  args[3] #CONFIG$vcfFileCommon
outPrefix <-  args[4] #"~/mutationSignatureGWAS"

mutationCols <- c("Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "fc")
covariateCols <- c("AGE", "GENDER", "BMI", "C1", "C2", "C3")
permute <- T
do_rankit_normalization <- T
yieldSize <- 2e6
genotypeMode <- "GT"
minorAlleleFreq <- 0.1
topNsignificant <- 100
ignoreCovariates <- c("ETHNCTY")

# For rSubmitter
GWASsource <- "~/scripts/FraserLab/somaticMutationsProject/R/GWAS.R"
extraBashLines <- c("module load mariadb/10.2.11", "module load R/3.4.0")
packages <- "VariantAnnotation"
workDir <- file.path("/scratch/users/paedugar/somaticMutationsProject/tempDir", "GWAS", paste0(sample(c(LETTERS, 0:9, 20, replace = T)), collapse = ""))
dir.create(workDir, recursive = T)

#permutations <- 100
#source("~/scripts/FraserLab/somaticMutationsProject/R/GWAS.R")
#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates_skin_fc/skin_C_T_fc_covariates.txt"
#library("rjson")
#CONFIG <- fromJSON(file = "../../../config.json")
#vcfFile <-  CONFIG$vcfFileCommon

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


#---------------------------
# MAIN

#####
# Reading mutation data
# Reading table
mutationMat_withCov <- read.table(mutationsFile, header = T, sep = "\t", stringsAsFactors = T)

# Select only caucasians
mutationMat_withCov <- mutationMat_withCov[mutationMat_withCov$RACE == 3,]

rownames(mutationMat_withCov) <- mutationMat_withCov$gtexId


# Select mutation types only
mutationMat <- t(mutationMat_withCov[,mutationCols])

if(do_rankit_normalization)
    mutationMat <- rankitNormalize(mutationMat)


#Read metadata
covariates <- t(mutationMat_withCov[,covariateCols])
covariates <- covariates[ !rownames(covariates) %in% ignoreCovariates, , drop = F] 


## Creating bins of genome
genomeRanges <- binGenome(CONFIG$auxiliaryFiles$hg19_chromSizes, yieldSize, exclude = c("Y", "X"))

################
# Goes over all mutation signatures

allSignatureGWAS <- list()
for (i in 1:nrow(mutationMat)) {
    phenoVector <- mutationMat[i,]

    # Makes permuted table 
    permuted <- do.call(rbind, lapply(rep(T, permutations), function(x,y) sample(y), y = phenoVector))
    # Real data is first row
    permuted <- rbind(phenoVector, permuted)

    # Doing GWAS
    #permutedResults <- lapply(genomeRanges[1:3], lapplyHelper, vcfFile = vcfFile, phenotype = permuted, genotype = genotypeMode, topNsignificant = topNsignificant, covariateMatrix = covariates, minorAlleleFreq = minorAlleleFreq, permute = permute)
    permutedResults <- "failed"
    tries <- 1
    while(is.character(permutedResults)) {
        cat("Running", tries, " try\n")
        permutedResults <- tryCatch(rSubmitter::superApply(genomeRanges, lapplyHelper, vcfFile = vcfFile, phenotype = permuted, 
                                      genotype = genotypeMode, topNsignificant = topNsignificant, 
                                      covariateMatrix = covariates, minorAlleleFreq = minorAlleleFreq, permute = permute,
                                      tasks = 200, mem = "4G", partition = "hbfraser,owners,hns", time = "1:30:00",
                                      sources = GWASsource, packages = packages, extraBashLines = extraBashLines, 
                                      workingDir =  workDir),
                                    error = function(e) { print(e); "failed" }
                                    )
        if (tries >= 4)
            stop("Reached 4 tries stopping now")
        tries <- tries + 1
    }

    mergedResults <- mergeGWASresults(permutedResults, topNsignificant)
    
    mergedResults$realTop$pvalue_perm <- c(sum(mergedResults$permutedMin[-1] <= mergedResults$realTop$pvalue[1] ) / length(mergedResults$permutedMin[-1]), rep(NA, topNsignificant-1))
    mergedResults$realTop$pvalue_bonf <- mergedResults$realTop$pvalue * mergedResults$totalSNPs
    mergedResults$realTop$signature <- rownames(mutationMat)[i]
    mergedResults$realTop$nTests <- mergedResults$totalSNPs
    
    allSignatureGWAS[[i]] <- mergedResults
}


pvals <- do.call(rbind, lapply(allSignatureGWAS, function(x) x$realTop))

# Writing results
write.table(pvals, paste0(outPrefix), sep = "\t", row.names = F, quote = F)
#save(allSignatureGWAS, file = "~/allSignatureGWAS.RData")

# Perfoms linear regression between given snps and mutation signatures from a H matrix
# do sh$ ml mariadb

library("MatrixEQTL")
library("ggplot2")
source("../../R/gtex.R", chdir = T)
source("../../R/GWAS.R", chdir = T)
source("../../R/misc.R", chdir = T)

cmdArgs <- commandArgs(T)

mutationSignature_H_file <- cmdArgs[1]
snps <- parseArg(cmdArgs[2], sep = "\t", trim = T)
vcfFile <- cmdArgs[3]

fromConsensus <- T

metadataCols <- c("SMRIN", #RIN
                  "SMTSISCH", #Ischemic time
                  "SMTSPAX" #Time in fixative
                  )

indInfo <- c("AGE", "GENDER", "BMI")
PCS <- c("C1", "C2", "C3")
#indInfo <- c("AGE", "ETHNCTY", "GENDER")

#mutationSignature_H_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Colon_Sigmoid-n6_0.0_0.7_H_consensus.txt"
#snps <- c("3_36395938_CA_C_b37", "3_36364118_C_A_b37")
#vcfFile <- "/scratch/PI/hbfraser/gtex/raw/Genotype/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz"

#-------------------
# METHODS

lmMatFunction <- function(x, y, useModel = modelLINEAR, cvrt = NULL, pvalCutoff = 1, errorCovariance = numeric(), outFile = tempfile(), min.pv.by.genesnp = F, noFDRsaveMemory = F){
    
    xMat <- SlicedData$new()
    yMat <- SlicedData$new()
    cvrtMat <- SlicedData$new()
    
    xMat$CreateFromMatrix(x)
    yMat$CreateFromMatrix(y)
    
    if(!is.null(cvrt)) 
        cvrtMat$CreateFromMatrix(cvrt)
    
    
    results <- Matrix_eQTL_engine(snps = xMat, gene = yMat, cvrt = cvrtMat, output_file_name = outFile,
                      pvOutputThreshold = pvalCutoff, useModel = useModel, min.pv.by.genesnp = min.pv.by.genesnp, noFDRsaveMemory = noFDRsaveMemory)
    
    return(results)
}

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

read_H_mat <- function(x, fromConsensus = F) {
    x <- read.table(x, stringsAsFactors = F, header = T, sep = "\t")
    if(fromConsensus) {
        rownames(x) <- x[,1]
        x <- x[,-1]
    } else {
        rownames(x) <- 1:nrow(x)
    }
    return(as.matrix(x))
}

#-------------------
# MAIN

Hmat <- read_H_mat(mutationSignature_H_file, fromConsensus = fromConsensus)
sraIds <- colnames(Hmat)

sraIds <- colnames(Hmat)
gtexIds <- sraToGtex(sraIds, formatOut = "long")
gtexIdsShort <- sraToGtex(sraIds, formatOut = "short")

# Read snps
genotype <- getSNPsByID(ids = snps, vcfFile = vcfFile, individuals = gtexIdsShort, genotype = "GT", chromFormat = "ENSEMBL")
genotype <- t(genotype)
genotype <- rangeSNPs(GRanges("chr3", IRanges(start = 37034841, end = 37092337)),500e3, 500e3, vcfFile, gtexIdsShort)
genotype <- t(genotype)


# Read covariates
sampleMeta <- t(readSampleAnnotationGTEX(metadataCols))
sampleMeta_shortIds <- gsub("(GTEX-\\w+?)-.+", "\\1", colnames(sampleMeta))

indMeta <- t(readMetadataGTEX(indInfo))
genotypePCA <- readGenotypePCA()[PCS,]

# Getting ids for which we have info for all
sharedInds <- gtexIdsShort %in% colnames(genotype) & 
    gtexIds %in% colnames(sampleMeta) &
    gtexIdsShort %in% colnames(indMeta) &
    gtexIdsShort %in% colnames(genotypePCA) 



sraIds <- sraIds[sharedInds]
gtexIds <- gtexIds[sraIds]
gtexIdsShort <- gtexIdsShort[sraIds]

# Creating metadata table
covariateMat <- rbind(indMeta[,gtexIdsShort], sampleMeta[,gtexIds])
covariateMat <- covariateMat[rowSums(!is.na(covariateMat)) == ncol(covariateMat),]
genotype <- genotype[,gtexIdsShort]
Hmat <- Hmat[,sraIds]

        
# Performing tests
Hmat_norm <- rankitNormalize(Hmat, 1)
x <- lmMatFunction(as.matrix(genotype), as.matrix(Hmat_norm), cvrt = as.matrix(covariateMat))
x <- x$all$eqtls
x$pvalue_bonf <- p.adjust(x$pvalue)
    


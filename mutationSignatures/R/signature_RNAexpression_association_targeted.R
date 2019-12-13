# Performs linear regressions between the expression of all genes and
# and the mutation signature loading for a given cancer type

library("MatrixEQTL")
library("ggplot2")
source("../../R/plots.R")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)

cmdArgs <- commandArgs(T)

pvalueCutoff <- as.numeric(cmdArgs[1])
mutationSignature_H_file <- cmdArgs[2]
targetGenes <- parseArg(cmdArgs[3], sep = ",", trim = T)
output <- cmdArgs[4]

metadataCols <- c("SMRIN", #RIN
                  "SMTSISCH", #Ischemic time
                  "SMTSPAX" #Time in fixative
                  )

fromConsensus <- T

indInfo <- c("AGE", "GENDER")
PCS <- c("C1", "C2", "C3", "C4")
#indInfo <- c("AGE", "ETHNCTY", "GENDER")





#pvalueCutoff <- 0.05
#mutationSignature_H_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Esophagus_Gastroesophageal_Junction-n6_0.0_0.7_H.txt"
#targetGenes <- parseArg("ENSG00000128383,ENSG00000076242,ENSG00000122512,ENSG00000139618,ENSG00000012048,ENSG00000179750,ENSG00000170734,ENSG00000154767,ENSG00000134574,ENSG00000140398", sep = ",", trim = T)
#geneNames <- "APOBEC3A,MLH1,PMS1,BRCA2,BRCA1,APOBEC3B,POLH,XPC,DDB2,NEIL1"
#output <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/axpressionAssociation_targeted/Esophagus_Gastroesophageal_Junction/n6_0.0_0.7_all.txt"

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



#-------------------
# MAIN

#pvalueCutoff <- 0.05
#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/"
#outDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/expressionAssociation/"
#H_files <- list.files(workingDir)[grep("H.txt", list.files(workingDir))]
#maf <- "n6_0.0_0.7"
#
#for(H_file in H_files) {
#tissue <- gsub("-n6_0.0_0.7_H.txt", "", H_file)
#mutationSignature_H_file <- file.path(workingDir, H_file)
#targetGenes <- parseArg("ENSG00000128383,ENSG00000076242,ENSG00000122512,ENSG00000139618,ENSG00000012048", sep = ",", trim = T)
#output <- file.path(outDir, tissue, "n6_0.0_0.7_all.txt")
#dir.create(dirname(outputSignificant), recursive = T)



# Read H matrix
Hmat <- read_H_mat(mutationSignature_H_file, fromConsensus = fromConsensus)
sraIds <- colnames(Hmat)
gtexIds <- sraToGtex(sraIds, formatOut = "long")
gtexIdsShort <- sraToGtex(sraIds, formatOut = "short")

# Reads expression file
expressionMat <- readAllGtexExpression(gtexIds)
rownames(expressionMat) <- gsub("\\.\\d+", "",  expressionMat[,1])
# SELECT GENES of interest
expressionMat <- expressionMat[rownames(expressionMat) %in% targetGenes,]


# Read covariates
sampleMeta <- t(readSampleAnnotationGTEX(metadataCols))
sampleMeta_shortIds <- gsub("(GTEX-\\w+?)-.+", "\\1", colnames(sampleMeta))

indMeta <- t(readMetadataGTEX(indInfo))
genotypePCA <- readGenotypePCA()[PCS,]

# Getting ids for which we have info for all
sharedInds <- gtexIds %in% colnames(expressionMat) & 
    gtexIds %in% colnames(sampleMeta) &
    gtexIdsShort %in% colnames(indMeta) &
    gtexIdsShort %in% colnames(genotypePCA) 



sraIds <- sraIds[sharedInds]
gtexIds <- gtexIds[sraIds]
gtexIdsShort <- gtexIdsShort[sraIds]



# Creating metadata table
covariateMat <- rbind(indMeta[,gtexIdsShort], sampleMeta[,gtexIds])
expressionMat <- expressionMat[,gtexIds]
# Only genes expressed in >20% of samples
#expressionMat <- expressionMat[rowSums(expressionMat > 1) >= (ncol(expressionMat) * 0.2),]
Hmat <- Hmat[,sraIds]
#colnames(covariateMat) <- sraIds
#colnames(expressionMat) <- sraIds

covariateMat <- covariateMat[rowSums(!is.na(covariateMat)) == ncol(covariateMat),]


for(i in 1:nrow(Hmat)) {
    
    # Selecting individuals that have non-zero values for mutation signatures
    inds <-  Hmat[i,] != 0
    phenoVector <- Hmat[ , inds]
    phenoVector <- rankitNormalize(phenoVector)[i,,drop =F]
    
    # Performing tests
    x <- lmMatFunction(phenoVector, as.matrix(expressionMat[,inds]), cvrt = as.matrix(covariateMat[,inds]))
    x <- x$all$eqtls
    x$pvalue_bonf <- p.adjust(x$pvalue)
    
    # Creating plots
    # Getting signatures normalized by covariates
    
    phenoVectNorm <- resid(lm(t(phenoVector) ~ t(covariateMat[,inds])))
    counter <- 1
    for(j in as.character(x$gene)) {
        if (counter > 10)
            break
        plot_title <- paste(c("Signature", rownames(Hmat)[i], counter, "Gene", j), collapse = "_")
        p <- scatter(data.frame(y = phenoVectNorm, x = as.numeric(expressionMat[j,inds])), x = "x", y = "y", regression =T)
        p <- p + ggtitle(plot_title) + xlab("TPM") + ylab("Normalized mutation signature weight")
        ggsave(file.path(dirname(output), paste0("plot_", plot_title, ".pdf")), p)
        counter <- counter + 1 

    }
    
    # Keeping results
    if(!exists("allHits")) {
        allHits <- x
    } else {
        allHits <- rbind(allHits, x)
    }
}


write.table(allHits, output, sep = "\t", col.names = T, row.names = F, quote = F)
rm(allHits)
#}

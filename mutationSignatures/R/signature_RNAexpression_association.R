# Performs linear regressions between the expression of all genes and
# and the mutation signature loading for a given cancer type

library("MatrixEQTL")
library("ggplot2")
source("../../R/plots.R")
source("../../R/gtex.R", chdir = T)
source("../../R/GO_analysis.R", chdir = T)
source("../../R/idConverter.R", chdir = T)

cmdArgs <- commandArgs(T)

pvalueCutoff <- as.numeric(cmdArgs[1])
mutationSignature_H_file <- cmdArgs[2]
outputSignificant <- cmdArgs[3]
outputNon <- cmdArgs[4]
outputGO <- cmdArgs[5]

fromConsensus <- T

metadataCols <- c("SMRIN", #RIN
                  "SMTSISCH", #Ischemic time
                  "SMTSPAX" #Time in fixative
                  )

indInfo <- c("AGE", "GENDER")
PCS <- c("C1", "C2", "C3")
#indInfo <- c("AGE", "ETHNCTY", "GENDER")





#pvalueCutoff <- 0.05
#mutationSignature_H_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Prostate-n6_0.0_0.7_H_consensus.txt"
#outputSignificant <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/expressionAssociation/Prostate/n6_0.0_0.7/significant.txt.gz"
#outputNon <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/expressionAssociation/Prostate/n6_0.0_0.7/nonSignificant.txt.gz"
#outputPlots <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/expressionAssociation/Prostate/n6_0.0_0.7/significant.pdf"
#outputGO <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/expressionAssociation/Prostate/n6_0.0_0.7/GO_results.txt"

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
#outputSignificant <- file.path(outDir, tissue, "n6_0.0_0.7_significant.txt.gz")
#outputNon <- file.path(outDir, tissue, "n6_0.0_0.7_nonSignificant.txt.gz")
#outputPlots <- file.path(outDir, tissue, "n6_0.0_0.7_significant.pdf")
#outputGO <- file.path(outDir, tissue, "GO_results.txt")
#dir.create(dirname(outputSignificant), recursive = T)



# Read H matrix
Hmat <- read_H_mat(mutationSignature_H_file, fromConsensus = fromConsensus)
sraIds <- colnames(Hmat)
gtexIds <- sraToGtex(sraIds, formatOut = "long")
gtexIdsShort <- sraToGtex(sraIds, formatOut = "short")

# Reads expression file
expressionMat <- readAllGtexExpression(gtexIds)
rownames(expressionMat) <- expressionMat[,1]

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
expressionMat <- expressionMat[rowSums(expressionMat > 1) >= (ncol(expressionMat) * 0.2),]
Hmat <- Hmat[,sraIds]
#colnames(covariateMat) <- sraIds
#colnames(expressionMat) <- sraIds

# Eliminate covariates missing or not variable
covariateMat <- covariateMat[rowSums(!is.na(covariateMat)) == ncol(covariateMat),]
mostRepeated <- apply(covariateMat,1,function(x) max(table(x)))
covariateMat <- covariateMat[mostRepeated/ncol(covariateMat) < 0.95,]

all_GO <- list()
all_GO_genes <- list()

for(i in 1:nrow(Hmat)) {
    
    # Selecting individuals that have non-zero values for mutation signatures
    inds <-  Hmat[i,] != 0
    phenoVector <- Hmat[ , inds]
    phenoVector <- rankitNormalize(phenoVector)[i,,drop =F]
    
    # Performing tests
    x <- lmMatFunction(phenoVector, as.matrix(expressionMat[,inds]), cvrt = as.matrix(covariateMat[,inds]))
    x <- x$all$eqtls
    x$pvalue_bonf <- p.adjust(x$pvalue)
    x$gene <- gsub("\\.\\d+", "", x$gene)
    x$alias <- ensemblToAlias_offline(x$gene)
    
    # Keeping results
    if(!exists("allHits")) {
        allHits <- x
    } else {
        allHits <- rbind(allHits, x)
    }
    
    # Skip if no significant results
    if (all(x$pvalue_bonf >= pvalueCutoff))
        next
    
    # Performing GO analyses
    GO_result <- performGO(setNames(x$pvalue_bonf, x$gene), gene_selection_function = function(p) p < pvalueCutoff)
    top_GO_categories <- get_GO_results(GO_result, fdr_signif=F, bonf_signif=F, n = 100)
    if(nrow(top_GO_categories) > 0) {
        top_GO_categories$signature <- rownames(Hmat)[i]
        
        # Get genes in to GO category
        genes_in_GO <- get_genes_GO(x$gene[x$pvalue_bonf < pvalueCutoff], top_GO_categories, GO_result, 15)
        genes_in_GO$alias <- ensemblToAlias_offline(genes_in_GO$gene)
        genes_in_GO$signature <-  rownames(Hmat)[i]
        
        all_GO[[i]] <- top_GO_categories
        all_GO_genes[[i]] <- genes_in_GO
    }
    
    # Creating plots
    # Getting signatures normalized by covariates
    
    phenoVectNorm <- resid(lm(t(phenoVector) ~ t(covariateMat[,inds])))
    counter <- 1
    significant <-  x[x$pvalue_bonf < pvalueCutoff,]
    for(k in 1:nrow(significant)) {
        j <- as.character(significant[k,"gene"])
        alias_gene <-  significant[k,"alias"]
        if (counter > 10)
            break
        plot_title <- paste(c("Signature", rownames(Hmat)[i], counter, "Gene", j, alias_gene), collapse = "_")
        p <- scatter(data.frame(y = phenoVectNorm, x = as.numeric(expressionMat[j,inds])), x = "x", y = "y", regression =T)
        p <- p + ggtitle(plot_title) + xlab("TPM") + ylab("Normalized mutation signature weight")
        ggsave(file.path(dirname(outputSignificant), paste0("plot_", plot_title, ".pdf")), p)
        counter <- counter + 1 

    }
    
}


all_GO <- do.call(rbind, all_GO)
all_GO_genes <- do.call(rbind, all_GO_genes)
write.table(all_GO_genes, paste0(outputGO, "_genes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(all_GO, outputGO, sep = "\t", col.names = T, row.names = F, quote = F)
write.table(allHits[allHits$pvalue_bonf < pvalueCutoff,], gzfile(outputSignificant), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(allHits[allHits$pvalue_bonf >= pvalueCutoff,], gzfile(outputNon), sep = "\t", col.names = T, row.names = F, quote = F)
rm(allHits)
#}

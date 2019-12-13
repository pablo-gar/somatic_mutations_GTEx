# Performs linear regressions between the expression of all genes and
# and mutation numbers

library("MatrixEQTL")
library("ggplot2")
source("../../R/plots.R")
source("../../R/gtex.R", chdir = T)
source("../../R/GO_analysis.R", chdir = T)
source("../../R/mutationGeneAnnotation.R", chdir = T)

cmdArgs <- commandArgs(T)

pvalueCutoff <- as.numeric(cmdArgs[1])
mutationsFile <- cmdArgs[2]
outputSignificant <- cmdArgs[3]
outputNon <- cmdArgs[4]
outputGO <- cmdArgs[5]

ignoreCovariates <- c("ETHNCTY")


#pvalueCutoff <- 0.05
#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount_include_distance_bias/mutation_covariates/Whole_Blood-n6_0.0_0.7.txt"
#outputSignificant <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutations_expressionAssociation/Whole_Blood/n6_0.0_0.7_significant.txt.gz"
#outputNon <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutations_expressionAssociation/Whole_Blood/n6_0.0_0.7_nonSignificant.txt.gz"
#outputPlots <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutations_expressionAssociation/Whole_Blood/n6_0.0_0.7_significant.pdf"
#outputGO <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutations_expressionAssociation/Whole_Blood/n6_0.0_0.7_GO_results.txt"

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
    
    rowNames <- rownames(x)
    colNames <- colnames(x)

    x <- apply(x, IND, function(x) qnorm((rank(x) - 0.5) / length(x)))
    if(IND == 1)
        x <- t(x)

    rownames(x) <- rowNames
    colnames(x) <- colNames
    return(x)

}

#-------------------



#-------------------
# MAIN

#####
# Reading mutation data
# Reading table
mutationMat_withCov <- read.table(mutationsFile, header = T, sep = "\t", stringsAsFactors = F)
mutationMat_withCov <- mutationMat_withCov[,!colnames(mutationMat_withCov) %in% c("C_C", "T_T")]

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
covariates <- covariates[ , colSums(is.na(covariates)) == 0 , drop = F] 

# Reads expression file
expressionMat <- readAllGtexExpression(gtexIdsLong)
rownames(expressionMat) <- expressionMat[,1]

# Getting ids for which we have info for all
sharedInds <- gtexIdsLong %in% colnames(expressionMat) & 
    gtexIdsLong %in% colnames(covariates) 

sraIds <- sraIds[sharedInds]
gtexIdsLong <- gtexIdsLong[sraIds]

# Only genes expressed in >20% of samples
expressionMat <- expressionMat[,gtexIdsLong]
expressionMat <- expressionMat[rowSums(expressionMat > 1) >= (ncol(expressionMat) * 0.2),]
expressionMat <- as.data.frame(rankitNormalize(as.matrix(expressionMat)))

mutationMat <- mutationMat[,gtexIdsLong]
covariates <- covariates[,gtexIdsLong]

# Eliminate covariates missing or not variable
covariates <- covariates[rowSums(!is.na(covariates)) == ncol(covariates),]
mostRepeated <- apply(covariates,1,function(x) max(table(x)))
covariates <- covariates[mostRepeated/ncol(covariates) < 0.95,]

all_GO <- list()
all_GO_genes <- list()

for(i in 1:nrow(mutationMat)) {
    
    # Selecting individuals that have non-zero values for mutation signatures
    inds <-  mutationMat[i,] != 0
    phenoVector <- mutationMat[ , inds]
    phenoVector <- rankitNormalize(phenoVector)[i,,drop =F]
    
    # Performing tests
    x <- lmMatFunction(phenoVector, as.matrix(expressionMat[,inds]), cvrt = as.matrix(covariates[,inds]))
    x <- x$all$eqtls
    x$pvalue_bonf <- p.adjust(x$pvalue)
    x$gene <- gsub("\\.\\d+", "", x$gene)
    x$alias <- queryRefCDS(x$gene, reference = "gene_id", query = "gene_name")
    
    
    # Keeping results
    if(!exists("allHits")) {
        allHits <- x
    } else {
        allHits <- rbind(allHits, x)
    }
    
    # Skip if no significant results
    if (all(x$pvalue_bonf >= pvalueCutoff))
        next
    
    ###
    # Performing GO analyses
    GO_result <- performGO(setNames(x$pvalue_bonf, x$gene), gene_selection_function = function(p) p < pvalueCutoff)
    top_GO_categories <- get_GO_results(GO_result, fdr_signif=F, bonf_signif=F, n = 100)
    
    if(nrow(top_GO_categories) > 0) {
        
        toPrintGO_n <- ifelse(nrow(top_GO_categories) > 15, 15, nrow(top_GO_categories))
        
        top_GO_categories$mut_type <- rownames(mutationMat)[i]
        
        # Get genes in to GO category
        genes_in_GO <- get_genes_GO(x$gene[x$pvalue_bonf < pvalueCutoff], top_GO_categories, GO_result, toPrintGO_n)
        genes_in_GO$alias <- queryRefCDS(genes_in_GO$gene, reference = "gene_id", query = "gene_name")
        genes_in_GO$mut_type <-  rownames(mutationMat)[i]
        
        all_GO[[i]] <- top_GO_categories
        all_GO_genes[[i]] <- genes_in_GO
    }
    
    ####
    # Creating plots association plots
    phenoVectNorm <- resid(lm(t(phenoVector) ~ t(covariates[,inds])))
    counter <- 1
    significant <-  x[x$pvalue_bonf < pvalueCutoff,]
    for(k in 1:nrow(significant)) {
        j <- as.character(significant[k,"gene"])
        alias_gene <-  significant[k,"alias"]
        if (counter > 10)
            break
        plot_title <- paste(c("mut_type", rownames(mutationMat)[i], counter, "Gene", j, alias_gene), collapse = "_")
        p <- scatter(data.frame(y = phenoVectNorm, x = as.numeric(expressionMat[j,inds])), x = "x", y = "y", regression =T)
        p <- p + ggtitle(plot_title) + xlab("TPM") + ylab("Normalized mutation load")
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

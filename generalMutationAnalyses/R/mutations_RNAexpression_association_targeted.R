# Performs linear regressions between the expression of all genes and
# and the mutation signature loading for a given cancer type

library("MatrixEQTL")
library("dplyr")
library("tidyr")
library("ggplot2")
library("metap")
source("../../R/FDR.R")
source("../../R/plots.R")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/mutationGeneAnnotation.R", chdir = T)

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

read_H_mat <- function(x) {
    x <- read.table(x, stringsAsFactors = F, header = T, sep = "\t")
    rownames(x) <- 1:nrow(x)
    return(as.matrix(x))
}

permuted_pvalues <- function(phenoVector, expressionMat, covariates, inds, targetGenes, nPerm) {
    
    # Helper function to do 1 permutation
    helper <- function(x, phenoVector, expressionMat, covariates, inds, targetGenes){
        targetGenes_random <- sample(rownames(expressionMat), length(targetGenes), replace = F)
        mat <- as.matrix(expressionMat[targetGenes_random,inds])
        rownames(mat) <- targetGenes
        result <- lmMatFunction(phenoVector, mat, cvrt = as.matrix(covariates[,inds])) 
        result <-  result$all$eqtls
        result$snps <- paste0(result$snps, ".", x)
        return(result)
    }
    
    # Permute
    eqtls <- lapply(1:nPerm, helper, phenoVector = phenoVector, expressionMat = expressionMat, covariates = covariates, inds = inds, targetGenes = targetGenes)
    eqtls <- do.call(rbind, eqtls)
    eqtls <- eqtls %>%
        dplyr::select(snps, gene, pvalue) %>%
        spread(gene, pvalue) 
    
    eqtls <- as.data.frame(eqtls)
    rownames(eqtls) <- eqtls$snps
    eqtls <- as.matrix(eqtls[,-1])
    
    return(eqtls)
    
}

# x is the results of linear regressions(data.frame)
# x_perm is results  from permuted_pvalues()
# returns fdr(or pvalue corrected by uniform expectation)
get_group_enrichment <- function(x, x_perm, targetGenes, pathway_gene_table) {
    
    real_pvalues <- setNames(x$pvalue, x$gene)
    
    results <- data.frame(group = unique(pathway_gene_table$pathway), pvalue = 1, stringsAsFactors = F)
    for(i in 1:nrow(results)) {
        group <- results$group[i]
        genes <- pathway_gene_table[pathway_gene_table$pathway == group, "gene"]
        genes <- genes[genes %in% colnames(x_perm)]
        
        perm_pvalues <- x_perm[,genes, drop = F]
        real <- real_pvalues[genes]
        
        if(ncol(perm_pvalues) == 1) {
            results[i, "pvalue"] <- sum(perm_pvalues[,1] <= real) / nrow(perm_pvalues)
        } else { 
            perm_better_than_real <- apply(perm_pvalues, 1, function(x) sumlog(x)$p) <= sumlog(real)$p
            results[i,"pvalue"]  <- sum(perm_better_than_real) / nrow(perm_pvalues)
        }
    }
    
    return(results)
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


all_group_enrichment <- list()
x_perm_all <- list()
for(i in 1:nrow(mutationMat)) {
    
    # Selecting individuals that have non-zero values for mutation signatures
    inds <-  mutationMat[i,] != 0
    phenoVector <- mutationMat[ , inds]
    phenoVector <- rankitNormalize(phenoVector)[i,,drop =F]
    
    # Performing tests
    #x <- lmMatFunction(phenoVector, as.matrix(expressionMat[targetGenes_ind,inds]), cvrt = as.matrix(covariates[,inds]))
    x <- lmMatFunction(phenoVector, as.matrix(expressionMat[,inds]), cvrt = as.matrix(covariates[,inds]))
    x_all <- x$all$eqtls
    x <- x_all[ x_all$gene %in% targetGenes, ]
    x$pvalue_bonf <- p.adjust(x$pvalue)
    x$alias <- queryRefCDS(x$gene, reference = "gene_id", query = "gene_name")
    x$pvalue_perm <- 1
    x$fdr <- FDRcalculation(x$pvalue, x_all$pvalue)
    
    # Performing permutations
    x_perm <- permuted_pvalues(phenoVector, expressionMat, covariates, inds, targetGenes, nPerm)
    x_perm_all[[ rownames(mutationMat)[i] ]] <- x_perm
    
    # Performs enrichment by group
    x_group_enrichment <- get_group_enrichment(x, x_perm, targetGenes, pathway_gene_table)
    x_group_enrichment$mut_type <- rownames(phenoVector)
    all_group_enrichment[[i]] <- x_group_enrichment
    
    
    ####
    # Creating plots association plots
    phenoVectNorm <- resid(lm(t(phenoVector) ~ t(covariates[,inds])))
    counter <- 1
    significant <-  x
    for(k in 1:nrow(significant)) {
        
        #Adding permutation pvalue
        gene_id <- as.character(significant[k, "gene"])
        x[k, "pvalue_perm"] <- sum(x[k, "pvalue"] >= x_perm[,gene_id]) / nPerm
        
        # Plotting scatter
        j <- as.character(significant[k,"gene"])
        alias_gene <-  significant[k,"alias"]
        plot_title <- paste(c("mut_type", rownames(mutationMat)[i], counter, "Gene", j, alias_gene), collapse = "_")
        p <- scatter(data.frame(y = phenoVectNorm, x = as.numeric(expressionMat[j,inds])), x = "x", y = "y", regression =T)
        p <- p + ggtitle(plot_title) + xlab("TPM") + ylab("Normalized mutation load")
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

all_group_enrichment <- do.call(rbind, all_group_enrichment)

write.table(all_group_enrichment, output_group_enrichment, sep = "\t", col.names = T, row.names = F, quote = F)
write.table(allHits, output, sep = "\t", col.names = T, row.names = F, quote = F)
save(x_perm_all, file = output_group_perm_tables)

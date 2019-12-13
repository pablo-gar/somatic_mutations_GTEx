
# and mutation numbers

source("../../R/ggthemes.R")
library("MatrixEQTL")
library("ggplot2")
library("dplyr")
library("tidyr")
source("association_functions.R")
source("../../R/plots.R")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/GO_analysis.R", chdir = T)
source("../../R/mutationGeneAnnotation.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    MUT_TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "all")
    
    cmdArgs <- commandArgs(T)

    mutation_map_file <- cmdArgs[1]
    metadata_file <- cmdArgs[2]
    tissue <- cmdArgs[3]
    phenotype_var <- cmdArgs[4]
    pvalueCutoff <- as.numeric(cmdArgs[5])
    outputSignificant <- cmdArgs[6]
    outputNon <- cmdArgs[7]
    outputGO <- cmdArgs[8]
    outputPlots <- cmdArgs[9]

    #pvalueCutoff <- 0.05
    #mutation_map_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/all_mutations_annotated.txt"
    #metadata_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/metadata.txt"
    #tissue <- "Liver"
    #phenotype_var <- "Nonsense"
    #outputSignificant <- file.path("/scratch/users/paedugar/somaticMutationsProject/vaf/expression_association_selection/", phenotype_var, tissue, "significant.txt.gz")
    #outputNon <- file.path("/scratch/users/paedugar/somaticMutationsProject/vaf/expression_association_selection", phenotype_var, tissue, "nonSignificant.txt.gz")
    #outputPlots <- file.path("/scratch/users/paedugar/somaticMutationsProject/vaf/expression_association_selection", phenotype_var, tissue, "plots")
    #outputGO <- file.path("/scratch/users/paedugar/somaticMutationsProject/vaf/expression_association_selection/", phenotype_var, tissue, "GO_results.txt")
    
    dir.create(outputPlots, recursive=T)
    # Read vaf and covariates
    mutation_map <- read.table(mutation_map_file, header=T, sep="\t", stringsAsFactors=F)
    mutation_map <- get_vaf_per_sample_impact(mutation_map, phenotype_var)
    mutation_map <- mutation_map[mutation_map$tissue == tissue,]

    covariates <- read_covariates(metadata_file, mutation_map$sraIds)
    
    # Selecting exactly the same individuals
    mutation_map <- mutation_map[mutation_map$sraIds %in% colnames(covariates),]
    covariates <- covariates[,mutation_map$sraIds]
    colnames(covariates) <- mutation_map$gtexIds_samples

    # Reads expression file
    expressionMat <- readAllGtexExpression(mutation_map$gtexIds_samples)
    rownames(expressionMat) <- gsub("\\.\\w+", "", expressionMat[,1])
    expressionMat <- expressionMat[,-1]

    # Getting ids for which we have info for all
    mutation_map <- mutation_map[mutation_map$gtexIds_samples %in% colnames(expressionMat),]
    covariates <- covariates[,mutation_map$gtexIds_samples]
    expressionMat <- expressionMat[,mutation_map$gtexIds_samples]

    # Only genes expressed in >20% of samples
    expressionMat <- expressionMat[rowSums(expressionMat > 1) >= (ncol(expressionMat) * 0.2),]
    expressionMat <- as.data.frame(rankitNormalize(as.matrix(expressionMat)))

    all_GO <- list()
    all_GO_genes <- list()

    for(i in MUT_TYPES[MUT_TYPES %in% colnames(mutation_map)]) {
        
        # Selecting individuals that have non-zero values for mutation signatures
        phenoVector <- t(as.matrix(setNames(mutation_map[,i,drop=T], mutation_map$gtexIds_samples)))
        phenoVector <- rankitNormalize(phenoVector, IND=1)
        
        # Performing tests
        x <- lmMatFunction(phenoVector, as.matrix(expressionMat), cvrt = as.matrix(covariates))
        x <- x$all$eqtls
        x$pvalue_bonf <- p.adjust(x$pvalue)
        x$gene <- gsub("\\.\\d+", "", x$gene)
        x$alias <- queryRefCDS(x$gene, reference = "gene_id", query = "gene_name")
        x$mut <- i
        x$tissue <- tissue 
        
        
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
            
            top_GO_categories$mut_type <- i
            
            # Get genes in to GO category
            genes_in_GO <- get_genes_GO(x$gene[x$pvalue_bonf < pvalueCutoff], top_GO_categories, GO_result, toPrintGO_n)
            genes_in_GO$alias <- queryRefCDS(genes_in_GO$gene, reference = "gene_id", query = "gene_name")
            genes_in_GO$mut_type <- i
            
            all_GO[[i]] <- top_GO_categories
            all_GO_genes[[i]] <- genes_in_GO
        }
        
        ####
        # Creating plots association plots
        phenoVectNorm <- resid(lm(t(phenoVector) ~ t(covariates)))
        counter <- 1
        significant <-  x[x$pvalue_bonf < pvalueCutoff,]
        for(k in 1:nrow(significant)) {
            j <- as.character(significant[k,"gene"])
            alias_gene <-  significant[k,"alias"]
            if (counter > 10)
                break
            plot_title <- paste(c("mut_type", i, counter, "Gene", j, alias_gene), collapse = "_")
            current_exp <- as.numeric(expressionMat[j,names(phenoVectNorm)])
            p <- scatter(data.frame(y = phenoVectNorm, x = current_exp), x = "x", y = "y", regression =T) + theme_sleek()
            p <- p + ggtitle(plot_title) + xlab("TPM") + ylab(paste0("VAF(Synonymous) / VAF(",phenotype_var, ")"))
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
    
}


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


main()

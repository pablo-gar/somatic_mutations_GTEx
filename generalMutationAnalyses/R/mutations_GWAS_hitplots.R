# Do module load mariadb/10.2.11
#
#' Takes output of GWAS script and makes boxplots of significant results

if(!dir.exists("../../R"))
    stop("Can't find the shared R path, make sure to run script from directory where it lives")

source("../../R/ggthemes.R", chdir = T)

library("ggplot2")
library("Homo.sapiens")
library("GenomicRanges")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")

source("../../R/gtex.R", chdir = T)
source("../../R/GWAS.R", chdir = T)
source("../../R/geneTools.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/plots.R", chdir = T)

GENE_ID_ALIAS <- as.data.frame(org.Hs.egSYMBOL)
GENE_ID_ENSEMBL <- as.data.frame(org.Hs.egENSEMBL)
GENE_ID_ENSEMBL <- GENE_ID_ENSEMBL[!duplicated(GENE_ID_ENSEMBL[,1]),]
rownames(GENE_ID_ALIAS) <- GENE_ID_ALIAS[,1]
rownames(GENE_ID_ENSEMBL) <- GENE_ID_ENSEMBL[,1]

# Gathering arguments
args <- commandArgs(T)
mutationsFile <- args[1] 
gwas_results_file <- args[2] 
vcfFile <-  args[3] #CONFIG$vcfFileCommon
extend_region <- as.numeric(args[4])
eqtl_file <- args[5]
outPrefix <-  args[6] #"~/mutationSignatureGWAS"

#outPrefix <-  "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_analyses/Brain_Hippocampus/n6_0.0_0.7_"
#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Brain_Hippocampus-n6_0.0_0.7.txt"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_rankit_caucasians/Brain_Hippocampus-n6_0.0_0.7.txt"
#outPrefix <-  "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_analyses/Brain_Nucleus_accumbens_basal_ganglia/n6_0.0_0.7_"

#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Brain_Nucleus_accumbens_basal_ganglia-n6_0.0_0.7.txt"
#eqtl_file <- "/scratch/PI/hbfraser/gtex/raw/Eqtls/eqtls/GTEx_Analysis_v7_eQTL_all_associations/Brain_Nucleus_accumbens_basal_ganglia.allpairs.txt.gz"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_rankit_caucasians/Brain_Nucleus_accumbens_basal_ganglia-n6_0.0_0.7.txt"

#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Brain_Hippocampus-n6_0.0_0.7.txt"
#eqtl_file <- "/scratch/PI/hbfraser/gtex/raw/Eqtls/eqtls/GTEx_Analysis_v7_eQTL_all_associations/Brain_Hippocampus.allpairs.txt.gz"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS/Brain_Hippocampus-n6_0.0_0.7.txt"
#outPrefix <-  "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_analyses/Brain_Hippocampus/n6_0.0_0.7_"

#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Breast_Mammary_Tissue-n6_0.0_0.7.txt"
#eqtl_file <- "/scratch/PI/hbfraser/gtex/raw/Eqtls/eqtls/GTEx_Analysis_v7_eQTL_all_associations/Breast_Mammary_Tissue.allpairs.txt.gz"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_rankit_caucasians/Breast_Mammary_Tissue-n6_0.0_0.7.txt"
#outPrefix <-  "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_rankit_caucasians/Breast_Mammary_Tissue-n6/n6_0.0_0.7_"

#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Adipose_Subcutaneous-n6_0.0_0.7.txt"
#eqtl_file <- "/scratch/PI/hbfraser/gtex/raw/Eqtls/eqtls/GTEx_Analysis_v7_eQTL_all_associations/Adipose_Subcutaneous.allpairs.txt.gz"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_no_rankit_caucasians/Adipose_Subcutaneous-n6_0.0_0.7.txt"

#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Artery_Tibial-n6_0.0_0.7.txt"
#eqtl_file <- "/scratch/PI/hbfraser/gtex/raw/Eqtls/eqtls/GTEx_Analysis_v7_eQTL_all_associations/Artery_Tibial.allpairs.txt.gz"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_no_rankit_caucasians/Artery_Tibial-n6_0.0_0.7.txt"

#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Skin_Not_Sun_Exposed_Suprapubic-n6_0.0_0.7.txt"
#eqtl_file <- "/scratch/PI/hbfraser/gtex/raw/Eqtls/eqtls/GTEx_Analysis_v7_eQTL_all_associations/Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_no_rankit_caucasians/Skin_Not_Sun_Exposed_Suprapubic-n6_0.0_0.7.txt"

#mutationsFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/mutation_covariates/Brain_Nucleus_accumbens_basal_ganglia-n6_0.0_0.7.txt"
#eqtl_file <- "/scratch/PI/hbfraser/gtex/raw/Eqtls/eqtls/GTEx_Analysis_v7_eQTL_all_associations/Brain_Nucleus_accumbens_basal_ganglia.allpairs.txt.gz"
#gwas_results_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS/Brain_Nucleus_accumbens_basal_ganglia-n6_0.0_0.7.txt"
#outPrefix <-  "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/GWAS_analyses/Brain_Nucleus_accumbens_basal_ganglia/n6_0.0_0.7_"

#extend_region <- 2e6
#library("rjson")
#CONFIG <- fromJSON(file = "../../config.json")
#vcfFile <-  CONFIG$vcfFileCommon


outFile <- paste0(outPrefix, "finished_plots.txt")
genotypeMode <- "GT"
minorAlleleFreq <- 0.1
ignoreCovariates <- c("ETHNCTY")
PVAL_SIGNIF <- 0.05

do_rankit_normalization <- T

#---------------------------
# METHODS

read_mutation_covariates <- function(mutationsFile, do_rankit_normalization) {
    
    mutationMat_withCov <- read.table(mutationsFile, header = T, sep = "\t", stringsAsFactors = T)
    sample_col <- which(colnames(mutationMat_withCov) ==  "sample")

    # Only include caucasian people
    mutationMat_withCov <- mutationMat_withCov[mutationMat_withCov[,"RACE"] == 3,]

    # Getting ids
    sraIds <- mutationMat_withCov$sample
    gtexIds <- sraToGtex(sraIds, formatOut = "short")
    gtexIds_long <- sraToGtex(sraIds, formatOut = "long")
    rownames(mutationMat_withCov) <- gtexIds


    # Select mutation types only
    mutationMat <- t(mutationMat_withCov[,1:(sample_col-1)])
    mutationMat <- mutationMat[!rownames(mutationMat) %in% c("C_C", "T_T"), ]
    
    if(do_rankit_normalization)
        mutationMat <- rankitNormalize(mutationMat)

    #Read metadata
    covariates <-t(mutationMat_withCov[(sample_col + 1) : ncol(mutationMat_withCov)])
    covariates <- covariates[ !rownames(covariates) %in% ignoreCovariates, , drop = F] 
    
    return(list(mutationMat, covariates, gtexIds_long))
        
}

convert_genotypes <- function(x) {
    
    key <- c(`0` = "Hom ref", `1` = "Het", `2` =  "Hom alt")
    x[!is.na(x)] <- key[as.character(x[!is.na(x)])]
    x <- factor(x, levels = key, ordered = T)
    return(x)
}

geno_id_to_string <- function(x, format_out = "string") {
    
    x <- strsplit(x, "_")
    if (format_out == "string") {
        x <- sapply(x, function (x) paste0("chr", x[1], ":", x[2], " / Ref: ", x[3], " / Alt: ", x[4]))
    } else if (format_out == "grange") {
        x <- lapply(x, function(x) data.frame(chr = paste0("chr", x[1]), start = as.numeric(x[2]), end = as.numeric(x[2]) + 1))
        x <- do.call(rbind, x)
    } else if (format_out == "grange_ensembl") {
        x <- lapply(x, function(x) data.frame(chr = x[1], start = as.numeric(x[2]), end = as.numeric(x[2]) + 1))
        x <- do.call(rbind, x)
    }
    return(x)
}

get_genes_around_snp <- function(snp, extend_region) {
    
    snp_pos <- makeGRangesFromDataFrame(geno_id_to_string(snp, format_out = "grange"))
    snp_range <- snp_pos + round(extend_region / 2)
    genes_around_snp <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), snp_range)
    genes_around_snp$gene_name <- GENE_ID_ALIAS[as.character(genes_around_snp$gene_id),2]
    genes_around_snp$gene_ensembl <- GENE_ID_ENSEMBL[as.character(genes_around_snp$gene_id),2]
    genes_around_snp$snp <- snp
    genes_around_snp$distance_tss <- 0
    genes_around_snp$description <- ensemblToDescription(genes_around_snp$gene_ensembl)
    
    genes_around_snp <- as.data.frame(genes_around_snp)
    
    genes_around_snp$distance_tss[genes_around_snp$strand == "-"] <- genes_around_snp$end[genes_around_snp$strand == "-"] - start(snp_pos)
    genes_around_snp$distance_tss[genes_around_snp$strand == "+"] <- genes_around_snp$start[genes_around_snp$strand == "+"] - start(snp_pos) 
    
    genes_around_snp <- genes_around_snp[order(genes_around_snp$distance_tss),]
    
    return(genes_around_snp)
}


plot_boxplot <- function(gwas_results, mutationMat, covariates, genotype_mat, mut_type, snp) {
    mutations <- resid(lm(mutationMat[,mut_type] ~ covariates))
    
    toPlot <- data.frame(mutations = mutations, genotype = convert_genotypes(genotype_mat[names(mutations),snp]))
    toPlot <- toPlot[!is.na(toPlot$genotype),]
    toPlotText <- paste0("   p = ", signif(gwas_results$pvalue, 2),
                         "; p_bonf = ", signif(gwas_results$pvalue_bonf, 2),
                         "\n   ", geno_id_to_string(snp),
                         "\n   r^2 = ", round((gwas_results$statistic / sqrt(gwas_results$degrees_freedom + gwas_results$statistic ^2) ) ^2, 2)
                         )
    ylabText <- paste0("Normalized mutations (", sub("_", ">", mut_type), ")")
    xlabText <- "Genotype"
    
    p <- ggplot(toPlot, aes(x = as.factor(genotype), y = mutations)) +
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(alpha = 0.3, width = 0.2) +
        annotate("text", label = toPlotText, x = -Inf, y = Inf, hjust = 0, vjust = 1, fontface = "italic") +
        xlab(xlabText) +
        ylab(ylabText) +
        theme_grid_y()
    
    return(p)
}


plot_expression_mut <- function(mutationMat, covariates, genotype_mat, mut_type, snp, expresion_mat, gene_id, gene_name) {
    
    mutations <- resid(lm(mutationMat[,mut_type] ~ covariates))
    sharedInd <- names(mutations)[names(mutations) %in% rownames(expression_mat)]
    
    if(length(sharedInd) == 0)
        return(NULL)
    if(!gene_id %in% colnames(expression_mat))
        return(NULL)
    
    toPlot <- data.frame(mutations = mutations[sharedInd], expression = expression_mat[sharedInd, gene_id])
    toPlot <- toPlot[toPlot$expression != 0 ,]
    toPlot$expression <- rankitNormalize_vector(toPlot$expression)
    #cor_results <- cor.test(toPlot[,1], toPlot[,2])
                        
                         
    ylabText <- paste0("Normalized mutations (", sub("_", ">", mut_type), ")")
    xlabText <- paste(gene_name, "expression normalized")
    
    p <- scatter(toPlot, x = "expression", y = "mutations", regression = T)
    p <- p + theme_noGrid() + xlab(xlabText) + ylab(ylabText)
    
    #toPlotText <- paste0("   p = ", signif(cor_results[["p.value"]], 2),
    #                     "\n   r^2 = ", round(cor_results[["estimate"]], 2)
    #                     )
    #p <- ggplot(toPlot, aes(x = as.factor(genotype), y = mutations)) +
    #    geom_point(outlier.shape = NA) + 
    #    geom_jitter(alpha = 0.3, width = 0.2) +
    #    annotate("text", label = toPlotText, x = -Inf, y = Inf, hjust = 0, vjust = 1, fontface = "italic") +
    #    xlab(xlabText) +
    #    ylab(ylabText) +
    #    theme_grid_y()
    
    return(p)
}


# Takes a gz file path and returns a read table after performing a zgrep
# on the file looking fo the snp,
find_eqtls <- function(eqtl_file, snp) {
    
    pattern <- paste(paste0("-e ", snp), collapse = " ")
    zgrep_results <- system(paste("zgrep", pattern, eqtl_file), intern = T)
    eqtls <- read.table(textConnection(zgrep_results), sep = "\t", stringsAsFactors = F, header = F)
    
    colnames(eqtls) <- c("gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se")
    
    eqtls <- eqtls[order(eqtls$pval_nominal), ]
    return(eqtls)
    
}

write_results_with_rsids <- function(gwas_results, out_prefix) {
    
    snp_granges <- geno_id_to_string(gwas_results$snps, format_out = "grange_ensembl")
    snp_rsid <- as.data.frame(snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, makeGRangesFromDataFrame(snp_granges)))
    ids <- paste0(snp_rsid$seqnames, ".", snp_rsid$pos)
    snp_rsid <- snp_rsid[!duplicated(ids),]
    rownames(snp_rsid) <- ids
    gwas_results$rsid <- snp_rsid[paste0(snp_granges$chr, ".", snp_granges$start),"RefSNP_id"]
    
    write.table(gwas_results, paste0(out_prefix, "gwas_results.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
}

plotEqtl <- function(gene_id, gene_name, snp, expression_mat, covariates, genotype_mat) {
    
    shared_ind <- rownames(genotype_mat) [rownames(genotype_mat) %in% rownames(expression_mat)]
    if(length(shared_ind) == 0)
        return(NULL)
        
    if(sum(colnames(expression_mat) %in% gene_id) != 1)
        return(NULL)
    
    # Selecting shared individuals
    genotype_mat <- genotype_mat[shared_ind,,drop=F]
    expression_mat <- expression_mat[shared_ind,,drop = F]
    
    # Transforming expression values and getting residuals after regression of covariates
    current_expression <- expression_mat[rownames(genotype_mat), gene_id, drop = F]
    current_expression <- current_expression[current_expression[,1] != 0, , drop = F]
    
    if(nrow(current_expression) / length(shared_ind) < 0.2)
        return(NULL)
    
    current_expression[,1] <- rankitNormalize_vector(current_expression[,1])
    
    aaa <- data.frame(covariates[rownames(current_expression), -(1:2)], expression = current_expression[,1])
    exp_resid <- tryCatch(resid(lm(expression ~., data = aaa)), error = function(e) NULL)
    
    if(is.null(exp_resid))
        return(NULL)
                              
    
    toPlot <- data.frame(genotype = factor(genotype_mat[names(exp_resid),snp]), expression = exp_resid)
    toPlot <- toPlot[!is.na(toPlot$genotype), ]
    
    corResults <- cor.test(as.numeric(toPlot$genotype), toPlot$expression)
    corResults <- data.frame(gene_id = gene_id, gene_name = gene_name, snp = snp, cor = corResults[["estimate"]], pvalue = corResults[["p.value"]])
    
    toPlot$genotype <- convert_genotypes(as.character(toPlot$genotype))
    toPlotText <- paste0("r^2 = ", round(corResults$cor[1] ^ 2, 2), "\np = ", signif(corResults$pvalue[1], 2))
    p <- ggplot(toPlot, aes(x = genotype, y = expression)) + 
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(alpha = 0.3, width = 0.2) +
        annotate("text", label = toPlotText, x = -Inf, y = Inf, hjust = 0, vjust = 1, fontface = "italic") +
        theme_grid_y()
    
    return(list(corResults, p))
    
}

#---------------------------
# MAIN

## Reading gwas results
gwas_results <- read.table(gwas_results_file, header = T, sep = "\t", stringsAsFactors = F)
gwas_results <- gwas_results[gwas_results$pvalue_bonf < PVAL_SIGNIF,]


# Stop if not signigicant results
if(nrow(gwas_results) == 0){
    writeLines("Done", outFile)
    quit("no")
}

# getting snp ranges
write_results_with_rsids(gwas_results, outPrefix)

# Get rsids


# Reading mutations and covariatres
mut_cov <- read_mutation_covariates(mutationsFile, do_rankit_normalization)
mutationMat <- mut_cov[[1]]
covariates  <- mut_cov[[2]]
gtexIds_long <- mut_cov[[3]]

# Reading genotypes
genotype_mat <- getSNPsByID(ids = gwas_results$snps, vcfFile = vcfFile, individuals = colnames(mutationMat), genotype = genotypeMode)

# Get only individuals for which we have data
mutationMat <- t(mutationMat[,rownames(genotype_mat)])
covariates <- t(covariates[,rownames(genotype_mat)])

# Read expression
expression_mat <- readAllGtexExpression(gtexIds_long)
rownames(expression_mat) <- gsub("\\..+", "", expression_mat[,1])
expression_mat <- expression_mat[,-1]
colnames(expression_mat) <- gsub("(GTEX-\\w+?)-.+", "\\1", colnames(expression_mat))
expression_mat <- t(expression_mat)

#####
# Ploting boxplots
eqtls <- find_eqtls(eqtl_file, gwas_results$snps)
genes_around_snp <- list()
eqtl_all <- list()
expression_mut_all <- list()
for(i in 1:nrow(gwas_results)) {
    
    mut_type <- gwas_results$signature[i]
    snp <- gwas_results$snps[i]
    
    # Making boxplot; TODO make this into a function
    p <- plot_boxplot(gwas_results[i,], mutationMat, covariates, genotype_mat, mut_type, snp)
    
    ggsave(paste0(outPrefix, "boxplot_", mut_type, "_", snp, ".pdf"), p, width = 3)
    
    # Finding closest genes
    genes_around_snp[[i]] <- get_genes_around_snp(snp, extend_region)
    
    # Doing own eqtl analyses
    for(j in 1:nrow(genes_around_snp[[i]])) {
        gene_id <- genes_around_snp[[i]]$gene_ensembl[j]
        gene_name <- genes_around_snp[[i]]$gene_name[j]
        id <- paste0(snp, gene_id)
        
        if(is.null(eqtl_all[[id]])) {
            eqtl_result <- plotEqtl(gene_id, gene_name, snp, expression_mat, covariates, genotype_mat)
            if(!is.null(eqtl_result)) {
                ggsave(paste0(outPrefix, "boxplot_eqtl", mut_type, "_", snp,  "_", gene_id, "_", gene_name, ".pdf"), eqtl_result[[2]], width = 3)
                eqtl_all[[id]] <- eqtl_result[[1]]
                
            }
        }
        
        id2 <- paste0(mut_type, gene_id)
        if(is.null(expression_mut_all[[id2]])) {
            # Expression association with mutaitons
            expresion_mut_plot <- plot_expression_mut(mutationMat, covariates, genotype_mat, mut_type, snp, expresion_mat, gene_id, gene_name)
            if(!is.null(expresion_mut_plot))
                ggsave(paste0(outPrefix, "scatter_expression_mutation_", mut_type, "_", gene_id, "_", gene_name, ".pdf"), expresion_mut_plot, height = 4, width = 4)
            expression_mut_all[[id2]] <- T
        }
    }
}
genes_around_snp <- do.call(rbind, genes_around_snp)

write.table(genes_around_snp, paste0(outPrefix, "genes_near_snps.txt"), sep = "\t", quote = F, col.names = T, row.names = T)
write.table(eqtls, paste0(outPrefix, "eqtls.txt"), sep = "\t", quote = F, col.names = T, row.names = T)
writeLines("Done", outFile)

if(length(eqtl_all) > 0) {
    eqtl_all <- do.call(rbind, eqtl_all)
    eqtl_all <- eqtl_all[order(eqtl_all$pvalue),]
    write.table(eqtl_all, paste0(outPrefix, "eqtls_mine.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
}

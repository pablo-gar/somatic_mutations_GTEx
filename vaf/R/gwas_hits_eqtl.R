# Do module load mariadb/10.2.11
#
#' Takes output of GWAS script and makes boxplots of significant results

suppressMessages({
    source("../../R/ggthemes.R", chdir = T)
    library("GenomicRanges")
    library("ggplot2")
    library("dplyr")
    library("tidyr")
    library("Homo.sapiens")
    library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    
    source("association_functions.R")
    source("../../R/gtex.R", chdir = T)
    source("../../R/GWAS.R", chdir = T)
    source("../../R/geneTools.R", chdir = T)
    source("../../R/misc.R", chdir = T)
    source("../../R/plots.R", chdir = T)
})

GENE_ID_ALIAS <- as.data.frame(org.Hs.egSYMBOL)
GENE_ID_ENSEMBL <- as.data.frame(org.Hs.egENSEMBL)
GENE_ID_ENSEMBL <- GENE_ID_ENSEMBL[!duplicated(GENE_ID_ENSEMBL[,1]),]
rownames(GENE_ID_ALIAS) <- GENE_ID_ALIAS[,1]
rownames(GENE_ID_ENSEMBL) <- GENE_ID_ENSEMBL[,1]

main <- function(cmdArgs = commandArgs(T)) {
    
    COVARIATES <- c('AGE', 'BMI', 'C1', 'C2', 'C3', 'RACE')
    MUT_TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "all")
    genotypeMode <- "GT"
    
    # Gathering arguments
    args <- commandArgs(T)
    snp_ids_file <- args[1] 
    mutation_map_file <- args[2]
    vcfFile <-  args[3] #CONFIG$vcfFileCommon
    phenotype_var <- args[4]
    tissue <- args[5] 
    extend_region <- as.numeric(args[6])
    outPrefix <-  args[7] #"~/mutationSignatureGWAS"

    #snp_ids_file <- '/scratch/users/paedugar/somaticMutationsProject/vaf/gwas_selection_eqtl/Nonsense/Artery_Aorta/signif_gwas_snps.txt'
    #mutation_map_file <- '/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/all_mutations_annotated.txt'
    #phenotype_var <- 'Nonsense'
    #vcfFile <- CONFIG$vcfFileCommon
    #tissue <- 'Artery_Aorta'
    #extend_region <- 2e6
    #covariate_file <- '/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/metadata.txt'
    #outPrefix <-  "/scratch/users/paedugar/somaticMutationsProject/vaf/gwas_clonality_eqtl/vaf_corrected/Muscle_Skeletal"


    #---------------------------
    # MAIN

    ## Reading mutation
    printTime('Reading mutation file \n')
    mutation_map <- read.table(mutation_map_file, header=T, sep="\t", stringsAsFactors=F)
    if(phenotype_var %in% c('Missense', 'Nonsense')) {
        mutation_map <- get_vaf_per_sample_impact(mutation_map, phenotype_var) 
    } else {
        mutation_map <- get_vaf_per_sample(mutation_map, phenotype_var) 
    }
    
    mutation_map <- mutation_map[mutation_map$tissue == tissue,]
    mutation_map$gtexIds <- gtexLongToShort(mutation_map$gtexIds_samples)
    
    printTime('DONE Reading mutation file \n')
    
    ## Reading gwas results
    snp_ids <- readLines(snp_ids_file)
    
    ## Get all individuals of a given tissues
    gtexIds_long <- tissueToGtexLong(tissue)[,2, drop=T]
    gtexIds <- unique(gtexLongToShort(gtexIds_long))

    # Reading genotypes
    genotype_mat <- getSNPsByID(ids = snp_ids, vcfFile = vcfFile, individuals = unique(gtexIds), genotype = genotypeMode)

    # Reading covariates
    sampleAnnotation <- readMetadataGTEX(c('AGE', 'BMI', 'GENDER'))
    genotypes <- t(readGenotypePCA()[1:3,])
    covariates <- transform(merge(sampleAnnotation, genotypes, by='row.names'), row.names = Row.names, Row.names = NULL)

    # Read expression
    printTime('Reading expression file \n')
    expression_mat <- readAllGtexExpression(gtexIds_long)
    rownames(expression_mat) <- gsub("\\..+", "", expression_mat[,1])
    expression_mat <- expression_mat[,-1]
    expression_mat <- t(expression_mat)
    
    printTime('DONE Reading expression file \n')
    #####
    # Doing eqtl analysis
    genes_around_snp <- list()
    eqtl_all <- list()
    mutation_expression_all <- list()
    expression_mut_all <- list()
    rsids <- gtexSnp_to_rsid(colnames(genotype_mat))
    
    for(i in 1:ncol(genotype_mat)) {
        
        snp <- colnames(genotype_mat)[i]
        rsid <- rsids[i]
        
        flush.console()
        printTime(paste0('Working with snp ', snp, '\n'))
        
        # Finding closest genes
        genes_around_snp[[i]] <- get_genes_around_snp(snp, extend_region)
        genes_around_snp[[i]] <- genes_around_snp[[i]][!is.na(genes_around_snp[[i]]$gene_ensembl),]
        
        # Doing own eqtl analyses
        for(j in 1:nrow(genes_around_snp[[i]])) {
            gene_id <- genes_around_snp[[i]]$gene_ensembl[j]
            gene_name <- genes_around_snp[[i]]$gene_name[j]
            id <- paste0(snp, gene_id)
            
            if(is.null(eqtl_all[[id]])) {
                eqtl_result <- plotEqtl(gene_id, gene_name, snp, expression_mat, covariates, genotype_mat)
                if(!is.null(eqtl_result)) {
                    # Save plot
                    ggsave(file.path(outPrefix, paste0("boxplot_eqtl", "_", snp, "_", rsid, "_", gene_id, "_", gene_name, ".pdf")), eqtl_result[[2]], width = 3)
                    eqtl_all[[id]] <- eqtl_result[[1]]
                    eqtl_all[[id]]$rsid <- rsid
                    
                    # Get expression correlation with mutation
                    mutation_expression_result <- plotMutations(gene_id, gene_name, expression_mat, mutation_map, MUT_TYPES, covariates)
                    
                    # Get genotype correlation with mutation
                    mutation_expression_result <- plotMutations(gene_id, gene_name, expression_mat, mutation_map, MUT_TYPES, covariates)
                    
                    if(!is.null(mutation_expression_result)) {
                        ggsave(file.path(outPrefix, paste0("scatter_mutation_expression", "_", snp, "_", rsid, "_", gene_id, "_", gene_name, ".pdf")), mutation_expression_result[[2]], width = 14, height=4)
                        mutation_expression_all[[id]] <- mutation_expression_result[[1]]
                    }
                }
            }
        }
    }
    
    genes_around_snp <- do.call(rbind, genes_around_snp)
    write.table(genes_around_snp, file.path(outPrefix, "genes_near_snps.txt"), sep = "\t", quote = F, col.names = T, row.names = T)
    
    if(length(eqtl_all) > 0) {
        eqtl_all <- do.call(rbind, eqtl_all)
        eqtl_all <- eqtl_all %>%
            group_by(gene_id) %>%
            mutate(p_bonf=p.adjust(pvalue), tissue=tissue)
        eqtl_all <- eqtl_all[order(eqtl_all$pvalue),]
        
        
        mutation_expression_all <- do.call(rbind, mutation_expression_all)
        mutation_expression_all$tissue <- tissue
        
        write.table(eqtl_all, file.path(outPrefix, "eqtls_mine.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        write.table(mutation_expression_all, file.path(outPrefix, "mutation_expression.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
    }
    
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

gtexSnp_to_rsid <- function(snp_ids) {
    
    snp_granges <- geno_id_to_string(snp_ids, format_out = "grange_ensembl")
    snp_rsid <- as.data.frame(snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, makeGRangesFromDataFrame(snp_granges)))
    ids <- paste0(snp_rsid$seqnames, ".", snp_rsid$pos)
    snp_rsid <- snp_rsid[!duplicated(ids),]
    rownames(snp_rsid) <- ids
    rsid <- snp_rsid[paste0(snp_granges$chr, ".", snp_granges$start),"RefSNP_id"]
    
    return(rsid)
    
}

plotEqtl <- function(gene_id, gene_name, snp, expression_mat, covariates = NULL, genotype_mat) {
    
    if(!gene_id %in% colnames(expression_mat))
        return (NULL)
    
    
    # Transforming expression values and getting residuals after regression of covariates
    current_expression <- process_expression(expression_mat, gene_id)
    
    # Getting genotype
    current_genotype <- as.data.frame(genotype_mat[,snp,drop=F])
    colnames(current_genotype) <- 'genotype'
    current_genotype$gtexIds <- rownames(current_genotype)
    
    # merging
    all_data <- left_join(current_expression, current_genotype, by='gtexIds')
    
    # adding covariates if available
    if(!is.null(covariates)) {
        covariate_cols <- colnames(covariates)
        covariates$gtexIds <- rownames(covariates)
        all_data <- left_join(all_data, covariates, by='gtexIds')
        all_data <- all_data[rowSums(is.na(all_data)) == 0,]
        
        all_data$exp_resid <- tryCatch(resid(lm(as.formula(paste0('expression ~ ', paste0(covariate_cols, collapse=' + '))), data=all_data)), 
                                       error = function(e) NULL
                                       )
    } else {
        all_data <- all_data[rowSums(is.na(all_data)) == 0,]
        all_data$exp_resid <- all_data$expression
    }
    
    if(is.null(all_data$exp_resid))
        return(NULL)
    
    
    toPlot <- data.frame(genotype = factor(all_data$genotype), expression = all_data$exp_resid)
    
    corResults <- tryCatch(cor.test(as.numeric(toPlot$genotype), toPlot$expression), error = function(e) NULL)
    
    if (!is.null(corResults)) {
        
        corResults <- data.frame(gene_id = gene_id, gene_name = gene_name, snp = snp, cor = corResults[["estimate"]], pvalue = corResults[["p.value"]])
        
        toPlot$genotype <- convert_genotypes(as.character(toPlot$genotype))
        toPlotText <- paste0("r = ", round(corResults$cor[1], 2), "\np = ", signif(corResults$pvalue[1], 2))
        p <- ggplot(toPlot, aes(x = genotype, y = expression)) + 
            geom_boxplot(outlier.shape = NA) + 
            geom_jitter(alpha = 0.3, width = 0.2) +
            annotate("text", label = toPlotText, x = -Inf, y = Inf, hjust = 0, vjust = 1, fontface = "italic") +
            theme_sleek() +
            theme_grid_y()
        
        return(list(corResults, p))
        
    } else {
        return(NULL)
    }
    
}

plotMutations <- function(gene_id, gene_name, expression_mat, mutation_map, mut_types, covariates) {
    
    if(!gene_id %in% colnames(expression_mat))
        return(NULL)
    
    
    # Transforming expression values and getting residuals after regression of covariates
    current_expression <- process_expression(expression_mat, gene_id)
    
    # Getting mutations
    mut_types <- mut_types[mut_types %in% colnames(mutation_map)]
    mutation_map[,mut_types] <- rankitNormalize(as.matrix(mutation_map[,mut_types]), 2)
    
    # merging with mutation map
    all_data <- left_join(current_expression, mutation_map, by='gtexIds')
    # adding covariates if available
    if(!is.null(covariates)) {
        covariate_cols <- colnames(covariates)
        covariates$gtexIds <- rownames(covariates)
        all_data <- left_join(all_data, covariates, by='gtexIds')
        all_data <- all_data[rowSums(is.na(all_data[,covariate_cols])) == 0,]
        
        all_data$exp_resid <- tryCatch(resid(lm(as.formula(paste0('expression ~ ', paste0(covariate_cols, collapse=' + '))), data=all_data)), 
                                       error = function(e) NULL
                                       )
    } else {
        all_data <- all_data[rowSums(is.na(all_data)) == 0,]
        all_data$exp_resid <- all_data$expression
    }
    
    if(is.null(all_data$exp_resid))
        return(NULL)
    
    # gathering metrics by mutation types
    all_data <- gather(all_data, key = 'mut_type', value = 'val', one_of(mut_types))
    
    
    corResults <- tryCatch(
                           {all_data %>%
                              group_by(mut_type) %>%
                              summarise(cor_mutation_metric = cor.test(val, expression)[['estimate']],
                                        pvalue_mutation_metric = cor.test(val, expression)[['p.value']],
                                        gene_name=gene_name,
                                        gene_id=gene_id) %>%
                              ungroup()},
                           error=function(err) NULL
                           )
    
    
    if(!is.null(corResults)) {
        
        p <- scatter(all_data, x='expression', y='val', facet_x='mut_type', 
                     alpha=1, scales='fixed', method_cor='pearson', regression=T) +
            theme_sleek()
            
        return(list(corResults, p))
        
    } else {
        return(NULL)
    }
    
    
}


# Selects one gene, does quantile normalization, renames column to 'expression'
# and replace long gtex names to short ones
process_expression <- function(expression_mat, gene_id) {
    
    current_expression <- as.data.frame(expression_mat[, gene_id, drop = F])
    colnames(current_expression) <- 'expression'
    current_expression <- current_expression[current_expression[,1] != 0, , drop = F]
    current_expression[,1] <- rankitNormalize_vector(current_expression[,1])
    current_expression$gtexIds <- gtexLongToShort(rownames(current_expression))
    
    return(current_expression)
}

main()

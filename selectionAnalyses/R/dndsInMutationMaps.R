#' Permforms dN/dS calculations on mutations
#' It takes > 1 mutation maps. It then performs dN/dS calculations
#' both global and local, genome-wide and gene-specific
#' It performs same calculations binning genes into different
#' average expression levels and binning mutations into different VAFs
#'
#' Positional parameter
#' @param mafBins - number of bins with equal size of mutations to create
#' @param minExpression - minimal average expression for a gene to be considered
#' @param quantileExpression - quantile used to divide highly vs lowly epxressed genes
#' @param outPrefix - prefix use for all output files
#' @param mutationMap1 [mutationMap2] [...] - mutation map files
#'
#' @return Creates a series of files that have this prefix:
#'     [outPrefix],geneExpression[TOP|BOTTOM|ALL_EXPRESSED],selectionMaf:[MAF_LOWER_BOUND]_[MAF_UPPER_BOUND]
#' And this suffices:
#'     sel_cv.txt - gene level dN/dS calculations using global estimators
#'     sel_loc.txt - gene level dN/dS calculations using local estimators
#'     globalInds.txt - genome-wide level dN/dS calculations using global estimators
#'     genesSynonymousMutations.txt - counts of number of synonymous mutations per gene
#'     genesMutations.txt - mutations per gene and their annotations
#'
#' @examples Must be run from folder where this script resides
#'     Rscript dndsInMutationMaps.R 3 1 0.8 liver /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.5
#'     Rscript dndsInMutationMaps.R 4 1 0.8 /scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Colon_Sigmoid/n6_0.0_0.5/dndsout_ /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Colon_Sigmoid/n6_0.0_0.5/*
#'
#'     Rscript dndsInMutationMaps.R 4 1 0.8 /scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood_EXO/n6_0.0_0.5/dndsout_ /scratch/users/paedugar/somaticMutationsProject/mutationCount/map_with_repeated_old/Whole_Blood_EXO/n6_0.0_0.7/*
    


library("reshape")
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
source("dndscvPablo.R")
source("../../R/mutationGeneAnnotation.R", chdir = T)
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)


#------------------------------------------
# PARS
#------------------------------------------

args <- commandArgs(TRUE)
mafBins <- as.numeric(args[1]) #3
minExpression <- as.numeric(args[2]) #1
quantileExpression <- as.numeric(args[3]) #0.8 # Gene with high expression are the ones that fall above this quantile
outPrefix <- args[4] #"skin"
mutationFiles <- args[5:length(args)]

doMedianExpression <- T
MIN_MUT <- 4
 
#mafBins <- as.numeric("4") #3
#minExpression <- as.numeric("1") #1
#quantileExpression <- as.numeric("0.8") #0.8 # Gene with high expression are the ones that fall above this quantile
#outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood/n6_0.0_0.5/dndsout_"
#
#tissueMaf <-"/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Whole_Blood/n6_0.0_0.7/"
#mutationFiles <- file.path(tissueMaf, list.files(tissueMaf))

#maf <- "n6_0.0_0.7"
#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map"
#mutationFiles <- vector()
#for(tissue in list.dirs(workingDir, recursive = F)[6]) {
#    tissueMaf <-  file.path(tissue, maf)
#    mutationFiles <- c(mutationFiles, file.path(tissueMaf, list.files(tissueMaf)))
#}


#------------------------------------------
# METHODS
#------------------------------------------

# Read all mutations from a tissue/maf
readMutationMaps <- function(mutationFiles) {
    
    mutationsAll <- list()
    for(i in 1:length(mutationFiles)) {
        
        mutationFile <- mutationFiles[i]
        
        currentSample <- gsub(".txt", "", basename(mutationFile))
        
        mutations <- read.table(mutationFile, sep = "\t", stringsAsFactors = F)
        
        colnames(mutations) <- c("chr", "pos", "ref", "mut", "context", "coverage", "nSupport")
        mutations$chr <- gsub("chr", "", mutations$chr)
        mutations$sample <- currentSample
        mutations$gtexId <- sraToGtex(currentSample, formatOut = "long")
        mutationsAll[[i]] <- mutations
        closeAllConnections()

    }

    mutationsAll <- do.call(rbind, mutationsAll)
    
    return(mutationsAll)
}

selecteGenesCoverage <- function(expressionTable, coverage, doMedian = T) {
    
    expressionTable[,"Name"] <- gsub("\\..+", "",  expressionTable[,"Name"])
    
    # MELT IDEA
    #expressionTable <- melt(expressionTable)
    #colnames(expressionTable) <- c("gene_id", "gtexId", "readCounts")
    #expressionTable$CDS_length <- queryRefCDS(expressionTable$gene_id, "gene_id", "CDS_length")
    #
    #expressionTable <- expressionTable[!is.na(expressionTable$CDS_length), ]
    ##expressionTable <- expressionTable[expressionTable$readCounts/expressionTable$CDS_length >= coverage,]
    #expressionTable <- expressionTable[expressionTable$readCounts > coverage,]
    #
    #
    #
    rownames(expressionTable) <- expressionTable[,1]
    expressionTable <- as.matrix(expressionTable[,-1])
    
    # Select only CDS genes
    expressionTable <- expressionTable[!is.na(queryRefCDS(rownames(expressionTable), "gene_id", "gene_name")),]
    
    if(doMedian){
        expressionTable <- apply(expressionTable,1,median)
        expressionTable <- matrix(expressionTable, ncol = 1, dimnames = list(names(expressionTable), "median"))
    }
         
    
    expressionTableLogical <- expressionTable > coverage
    # Gets all by all values
    genes <- list()
    inds <- list()
    counts <- list()
    for(col in 1:ncol(expressionTable)) {
        nGenes <- sum (expressionTableLogical[,col] > 0)
        if(nGenes > 0) {
            inds[[col]] <- rep(colnames(expressionTableLogical)[col], nGenes)
            genes[[col]] <- rownames(expressionTableLogical)[expressionTableLogical[,col]]
            counts[[col]] <- expressionTable[expressionTableLogical[,col], col]
        }
    }
    results <- data.frame(gene_id = unlist(genes), gtexId = unlist(inds), readCounts = unlist(counts), stringsAsFactors = F)
    
    
    
    return(results)
}

# Gets mean and confidence intervals from a dndns table of the indicate columns
get_dnds_mean <- function(x, min_mut, col_names) {
    x <- x[rowSums(x[,2:5]) >= min_mut,]
    
    all_results <- list()
    
    funs <- list(mean = function(x) mean (x, trim = 0.1), median = median)
    
    for(j in seq_along(funs)) {
        
        fun <- funs[[j]]
        fun_name <- names(funs)[j]
        
        results <- data.frame(name = paste0(col_names, "_", fun_name), mle = 0, cilow = 0, cihigh = 0)
        for(i in seq_along(col_names)) {
            
            confidence_intervals <-  bootstrap_confidence_interval(x[,col_names[i]], fun, bootstrap_size=1, bootstrap_counts=10000)
            results[i, "mle"] <- fun(x[,col_names[i]])
            results[i, "cilow"] <- confidence_intervals[1]
            results[i, "cihigh"] <- confidence_intervals[2]
        }
        all_results <- c(all_results, list(results))
    }
    
    all_results <- do.call(rbind, all_results)
    rownames(all_results) <- 1:nrow(all_results)
    
    return(all_results)
}
    

    
#------------------------------------------
# MAIN
#------------------------------------------

# Read mutations and annotate
mutationsAll <- readMutationMaps(mutationFiles)


# Get mutation MAF quantiles
MAF <- mutationsAll$nSupport/mutationsAll$coverage
quantMaf <- quantile(MAF,seq(0,1, length.out = mafBins + 1))



# Read gene counts
#geneReadCounts <- readAllGeneReadCounts(unique(mutations$gtexId))
#geneReadCounts <- selecteGenesCoverage(geneReadCounts, 40)


# Do not work with expression if working with EXOME DATA

if(any(grepl("EXO", mutationFiles))) {

    genes <- list()
    genes[["AllExpressed"]] <- queryRefCDS(referenceField = "gene_id", queryField = "gene_name")
    
} else {
    
    geneReadCounts <- readAllGtexExpression(unique(mutationsAll$gtexId))
    geneReadCountsExpressed <- selecteGenesCoverage(geneReadCounts, minExpression, doMedianExpression)
    quantGenes <- quantile(geneReadCountsExpressed$readCounts, c(quantileExpression))

    # Gets top expressed genes
    genes <- list()
    genes[["Top"]] <- queryRefCDS(unique(geneReadCountsExpressed[geneReadCountsExpressed$readCounts > quantGenes[1], "gene_id"]), "gene_id", "gene_name")
    genes[["Bottom"]] <- queryRefCDS(unique(geneReadCountsExpressed[geneReadCountsExpressed$readCounts < quantGenes[1], "gene_id"]), "gene_id", "gene_name")
    genes[["AllExpressed"]] <- queryRefCDS(unique(geneReadCountsExpressed$gene_id), "gene_id", "gene_name")
    
}


# Do selection analyses
#for(i in "AllExpressed") {
for(i in names(genes)) {
    for(j in 1:(length(quantMaf))) {
        
        if(j < length(quantMaf)) {
            MAFlower <- quantMaf[j]
            MAFupper <- quantMaf[j+1]
            MAFlowerLabel <- names(quantMaf)[j]
            MAFupperLabel <- names(quantMaf)[j+1]
        } else { #last iteration is across entire MAF spectrum
            MAFlower <- 0
            MAFupper <- 0.5
            MAFlowerLabel <- "0%"
            MAFupperLabel <- "100%"
        }
        mutations <- mutationsAll[MAF >= MAFlower & MAF < MAFupper,]
        
        if(nrow(mutations) == 0)
            next
        
        dndsout <- dndscv(mutations[,c("sample", "chr", "pos", "ref", "mut")], gene_list = genes[[i]], cv = NULL, max_coding_muts_per_sample = 40e3, max_muts_per_gene_per_sample = 1000)
        
        currentPrefix <- paste0(outPrefix,",geneExpression:", i , ",selectionMaf:", MAFlowerLabel, "_", MAFupperLabel) 
        
        dndsloc_mean <- get_dnds_mean(dndsout$sel_loc, min_mut = MIN_MUT, col_names = c("wmis_loc", "wnon_loc"))
        dndsout$globaldnds <- rbind(dndsout$globaldnds, dndsloc_mean)
        
        write.table(dndsout$sel_cv, paste0(currentPrefix, ",_sel_cv.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        write.table(dndsout$sel_loc, paste0(currentPrefix, ",_sel_loc.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        write.table(dndsout$globaldnds, paste0(currentPrefix, ",_globalInds.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        write.table(as.data.frame(sort(table(dndsout$annotmuts$gene[dndsout$annotmuts$impact=="Synonymous"]),decreasing=T)),
                    paste0(currentPrefix, ",_genesSynonymousMutations.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        write.table(dndsout$annotmuts, paste0(currentPrefix, ",_genesAllMutations.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
    }
}

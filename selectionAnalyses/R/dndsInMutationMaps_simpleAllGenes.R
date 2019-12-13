#' Permforms dN/dS calculations on mutations
#' It takes > 1 mutation maps, and a table with gene symbols in the first column 
#' It then performs dN/dS calculations on that set of genes
#' both global and local, genome-wide and gene-specif
#'
#' Positional parameter
#' @param outPrefix - prefix use for all output files
#' @param gene_file - path to tab-separated file where first column has gene symbols 
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


library("reshape")
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
source("dndscvPablo.R")
source("../../R/mutationGeneAnnotation.R", chdir = T)
source("../../R/gtex.R", chdir = T)


#------------------------------------------
# PARS
#------------------------------------------

args <- commandArgs(TRUE)
outPrefix <- args[1]
mutationFiles <- args[2:length(args)]

 
#outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/selection/dndsout_tissuesMerged/n6_0.0_0.7/dndsout"
#tissueMaf <-"/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/"
#mutationFiles <- file.path(tissueMaf, list.files(tissueMaf,recursive = T))


#------------------------------------------
# METHODS
#------------------------------------------

# Read all mutations from a tissue/maf
readMutationMaps <- function(mutationFiles) {
    
    mutationsAll <- list()
    for(i in 1:length(mutationFiles)) {
        
        flush.console()
        print(mutationFiles[i])
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

    
#------------------------------------------
# MAIN
#------------------------------------------

# Read mutations and annotate
mutationsAll <- readMutationMaps(mutationFiles)

# Get genes of interest that are presente in dndscv
dndsout <- dndscv(mutationsAll[,c("sample", "chr", "pos", "ref", "mut")], cv = NULL, max_coding_muts_per_sample = 40e3, max_muts_per_gene_per_sample = 1000)


write.table(dndsout$sel_cv,  paste0(outPrefix, "dndsout_sel_cv.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(dndsout$sel_loc,  paste0(outPrefix, "dndsout_sel_loc.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(dndsout$globaldnds,  paste0(outPrefix, "dndsout_globalInds.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

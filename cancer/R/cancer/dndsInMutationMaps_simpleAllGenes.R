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
source("../../../selectionAnalyses/R/dndscvPablo.R")
source("../../../R/mutationGeneAnnotation.R", chdir = T)
source("../../../R/gtex.R", chdir = T)
source("../../../R/misc.R", chdir = T)


#------------------------------------------
# PARS
#------------------------------------------

args <- commandArgs(TRUE)
outPrefix <- args[1]
mutations_to_eliminate_file <- args[2]
mutationFiles <- args[3:length(args)]

 
#outPrefix <- "/scratch/users/paedugar/somaticMutationsProject/cancer/dndsout_tissuesMerged_noFalsePos/n6_0.0_0.7/"
#tissueMaf <-"/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/"
#mutationFiles <- file.path(tissueMaf, list.files(tissueMaf,recursive = T))
#mutationFiles <- mutationFiles[!grepl("EXO", mutationFiles)]
#mutations_to_eliminate_file <- "/scratch/users/paedugar/somaticMutationsProject/cancer/driverMutations/mutations_false_positives/n6_0.0_0.7/cancer_mutations_false_positives.txt"


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

get_dnds_mean <- function(dnds, get_log2 = T, n_mutations_cutoff = 5) {

    # Getting wmin average and non sense
    dnds_signif_mis <- get_dnds_good_genes(dnds, "wmis_loc", n_mutations_cutoff = n_mutations_cutoff)$wmis_loc
    dnds_signif_non <- get_dnds_good_genes(dnds, "wnon_loc", n_mutations_cutoff = n_mutations_cutoff)$wnon_loc

    if(get_log2) {
        dnds_signif_mis <- log2(dnds_signif_mis)
        dnds_signif_non <- log2(dnds_signif_non)
    }   

    dnds_mis_conf <- bootstrap_confidence_interval(dnds_signif_mis, mean)
    dnds_non_conf <- bootstrap_confidence_interval(dnds_signif_non, mean)

    dnds_mis_conf_median <- bootstrap_confidence_interval(dnds_signif_mis, median)
    dnds_non_conf_median <- bootstrap_confidence_interval(dnds_signif_non, median)

    dnds_mis_mean <- data.frame(mean = mean(dnds_signif_mis), conf_lo = dnds_mis_conf[1], conf_up = dnds_mis_conf[2],
                                median = median(dnds_signif_mis), conf_lo_median = dnds_mis_conf_median[1], conf_up_median = dnds_mis_conf_median[2]
                                )
    dnds_non_mean <- data.frame(mean = mean(dnds_signif_non), conf_lo = dnds_non_conf[1], conf_up = dnds_non_conf[2],
                                median = median(dnds_signif_non), conf_lo_median = dnds_non_conf_median[1], conf_up_median = dnds_non_conf_median[2]
                                )

    return(list( mis = dnds_mis_mean, non = dnds_non_mean))

}

get_dnds_good_genes <- function(dnds, measure = "wmis_loc", n_mutations_cutoff = 5) {
    
    if(measure == "wmis_loc") {
        mutColumn <- "n_mis"
        #dnds[,measure] <- dnds[p.adjust(dnds$pmis_loc, method = "BH") < 0.01, measure] 
    } else {
        mutColumn <- "n_non"
    }

    # only mutations with at least one syn and one from measure of interest
    dnds <- dnds[dnds$n_syn > 0 & dnds[,mutColumn] > 0, ]

    # Only genes that have the number of mutations cutoff
    dnds <- dnds[dnds$n_syn +  dnds[,mutColumn] >= n_mutations_cutoff, ]

    # if measure equals 0 lets create the lowest available
    dnds[dnds[,measure] == 0, measure] <- 1 / (max(dnds[,measure]) + 1)

    return(dnds)

}

    
#------------------------------------------
# MAIN
#------------------------------------------

# Read mutations and annotate
mutationsAll <- readMutationMaps(mutationFiles)
backup <- mutationsAll

# Eliminate mutations that are not desirible (i.e. useful for cancer analyses)
mutations_to_eliminate <- read.table(mutations_to_eliminate_file, sep = "\t", stringsAsFactors = F, header = F)
if(ncol(mutations_to_eliminate) > 4)
    mutations_to_eliminate <- mutations_to_eliminate[,1:4]

# Eliminate mutations
all_mut_ids <- paste0(mutationsAll[,1], ".", mutationsAll[,2], ".", mutationsAll[,3], ".", mutationsAll[,4])
mut_to_eliminate <- rep(F, nrow(mutationsAll))
for(i in 1:nrow(mutations_to_eliminate)) {
    mut_id <- paste(mutations_to_eliminate[i,], collapse = ".")
    is_mut_id <- all_mut_ids == mut_id
    flush.console()
    print(sum(is_mut_id))
    mut_to_eliminate <- mut_to_eliminate | is_mut_id
}

mutationsAll <- mutationsAll[!mut_to_eliminate,]


# Get genes of interest that are presente in dndscv
dndsout <- dndscv(mutationsAll[,c("sample", "chr", "pos", "ref", "mut")], cv = NULL, max_coding_muts_per_sample = 40e3, max_muts_per_gene_per_sample = 1000)

median_dnds <- get_dnds_mean(dndsout$sel_loc)


write.table(dndsout$sel_cv,  paste0(outPrefix, "dndsout_sel_cv.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(dndsout$sel_loc,  paste0(outPrefix, "dndsout_sel_loc.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(dndsout$globaldnds,  paste0(outPrefix, "dndsout_globalInds.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

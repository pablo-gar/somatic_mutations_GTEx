# Performs basic analysis of VAF and metadata
source("../../R/gtex.R", chdir=T)
source("../../R/plots.R", chdir=T)
source("../../R/misc.R", chdir=T)
library("dplyr")
library("tidyr")
library("ggplot2")

main <- function(cmdArgs = commandArgs(T)) {
    
    mutation_map_file <- cmdArgs[1]
    annotated_mut_file <- cmdArgs[2]
    cancer_mutation_dir <- cmdArgs[3]
    maf <- cmdArgs[4]
    out_file_map <- cmdArgs[5]
    metadata_out_file <- cmdArgs[6]
    
    #mutation_map_file <- "/scratch/users/paedugar/somaticMutationsProject/supp_info/tables/Table_mutations_all.tsv.gzip.temp.original"
    #annotated_mut_file <- "/scratch/users/paedugar/somaticMutationsProject/cancer/dndsout_tissuesMerged_all/1_all_mutations_annotated_dnds.txt.gz"
    #out_file_map <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/all_mutations_annotated.txt"
    #metadata_out_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/metadata.txt"
    #cancer_mutation_dir  <- "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/map/"
    #maf  <- "n6_0.0_0.7"
    
    ## READ mutations
    mutation_map <- read.table(mutation_map_file, sep="\t", stringsAsFactors=F, header=T)
    
    mutation_map$mut <- get_mutation(mutation_map)
    mutation_map$vaf <- with(mutation_map, alt_count/coverage)
    mutation_map$vaf_corrected <- rolling_median(mutation_map$vaf, mutation_map$coverage, windows=500, by_quantile=F)
    mutation_map$id <- paste0(mutation_map$sraIds, ".", mutation_map$chr, ".", mutation_map$pos, ".", mutation_map$ref, ".", mutation_map$alt)
    
    ## READ annotated mutations
    annotated_mutations <- read.table(annotated_mut_file, header = T, sep = "\t", stringsAsFactors = F)
    annotated_mutations <- annotated_mutations[,c("chr", "pos", "ref", "mut", "impact", "gene_id", "gene_name", "cds_length", "strand", "aachange", "ntchange", "sample")]
    annotated_mutations <- unique(annotated_mutations)
    annotated_mutations$id <- paste0(annotated_mutations$sample, ".chr", annotated_mutations$chr, ".", annotated_mutations$pos, ".", annotated_mutations$ref, ".", annotated_mutations$mut)
    annotated_mutations <- annotated_mutations[,c("id", "impact", "gene_id", "gene_name", "cds_length", "strand", "aachange", "ntchange")]
    
    mutation_map <- left_join(mutation_map, annotated_mutations, by="id")
    
    # READ cancer muattions
    mutations_cancer <- read_cosmic_maps(cancer_mutation_dir, maf)
    mutations_cancer$id <- paste0(mutations_cancer$sample, ".", mutations_cancer$chr, ".", mutations_cancer$start, ".", mutations_cancer$ref, ".", mutations_cancer$alt)
    mutations_cancer <- mutations_cancer[,c("id", colnames(mutations_cancer)[grep("cosmic", colnames(mutations_cancer))])]
    
    mutation_map <- left_join(mutation_map, mutations_cancer, by = "id")
    mutation_map$cancerous <- ifelse(mutation_map$alt_cosmic == mutation_map$alt, "Seen_in_cancer", "Not_in_cancer")
    
    ## READ metadata
    metadata <- appendMetadata(mutation_map)
    
    ## WRITE
    write.table(mutation_map, out_file_map, sep="\t", quote=F, row.names=F)
    write.table(metadata, metadata_out_file, sep="\t", quote=F, row.names=F)
    
    
}

appendMetadata <- function(x) {
    
    readCountFile = CONFIG[["auxiliaryFiles"]][["readCountTable"]]
    transDiversityFile = file.path(CONFIG[["auxiliaryFiles"]][["dir"]], "transcriptome_shannon_diversity.txt")
    
    metadata <-  read_all_metadata(unique(x$sraIds), readCountFile, transDiversityFile, indInfo = c("AGE", "GENDER", "BMI", "RACE"),  na.rm = F)
    metadata$sraIds <- rownames(metadata)
    
    return(metadata)
    
}


read_cosmic_maps <- function (mutationDir, maf) {

    results <- list()
    i <- 1
    for(tissue in list.dirs(mutationDir, recursive = F)) {
        if (grepl("EXO", tissue))
            next
        for(current_sample in list.files(file.path(tissue, maf))) {
            info <- file.info(file.path(tissue, maf, current_sample))
            if(any(info$size == 0))
                next

            current <- read.table(file.path(tissue, maf, current_sample), sep = "\t", stringsAsFactors = F, header = F)
            colnames(current) <- c("chr_cosmic", "start_cosmic", "end_cosmic", "ref_cosmic", "alt_cosmic", "strand", "chr", "start",
                                   "end", "ref", "alt", "context", "coverage", "alt_count")
            current$tissue <- basename(tissue)
            current$sample <- gsub(".txt", "", basename(current_sample))
            results[[i]] <- current
            i <- i + 1
        }
    }

    results <- do.call(rbind, results)

    return(results)

}

get_mutation <- function(x) {
    
    nuc <- c(A="T", T="A", C="G", G="C")
    
    from <- x$ref
    to <- x$alt
    
    to[from %in% c("A", "G")] <- nuc[to[from %in% c("A", "G")]]
    from[from %in% c("A", "G")] <- nuc[from[from %in% c("A", "G")]]
    
    mut <- paste0(from, ">", to)
    return(mut)
    
}

main()

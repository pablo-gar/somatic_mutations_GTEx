main <- function(cmdArgs = commandArgs(T)) {
    working_dir <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    x <- read_as_maf(working_dir)
    
    write.table(x, out_file, sep = "\t", quote = F, row.names = F, col.names = T)
}
                 
# Read mutations from a dir containing maps from the gtex project 
read_as_maf <- function(x) {
    
    sams <- list.files(x, full = T)
    
    results <- list()
    for(sam in sams) {
        current <- read.table(sam, sep = "\t", stringsAsFactors = F)
        colnames(current) <- c("chr", "pos", "Reference_Allele", "Tumor_Seq_Allele2", "ref_context", "coverage", "alt_counts")
        current$Variant_Type <- "SNP"
        current$Tumor_Sample_Barcode <- gsub(".txt", "", basename(sam))
        results[[sam]] <- current
    }
    
    results <- do.call(rbind, results)
    rownames(results) <- 1:nrow(results)
    
    return(results)
    
}

main()

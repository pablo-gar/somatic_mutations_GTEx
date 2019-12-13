main <- function(cmdArgs = commandArgs(T)) {
    
    output_prefix <- cmdArgs[1]
    mutation_table_files <- cmdArgs[-1]
    
    mutation_table <- read_tables(mutation_table_files)
    
    for(gene in unique(mutation_table$gene_name)) {
        
        current <- mutation_table[mutation_table$gene_name == gene,]
        current <- current[ , colnames(current) != "gene_name"]
        
        write.table(current, paste0(output_prefix, gene, ".txt"), sep = "\t", quote = F, row.names = F, col.names = T)
    }
    
}

read_tables <- function(x) {
    
    results <- list()
    
    for(i in x) {
        current <- read.table(i, sep = "\t", stringsAsFactors = F, header = T)
        tissue <- gsub(".txt", "", basename(i))
        current$tissue <- tissue
        current <- current[ , c("sample", "tissue", "chr", "pos", "pos", "ref", "mut", "gene_name")]
        results[[i]] <- current
        
    }
    
    results <- do.call(rbind, results)
    colnames(results) <- c("Sample_ID", "Cancer_Type", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Variant_Allele", "gene_name")
                           
    
    return(results)
    
}

main()

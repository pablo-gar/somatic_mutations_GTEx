# converts decode_DNMs.tsv to individual maps of mother and father mutations

main <- function(cmdArgs = commandArgs(T)){
    
    input <- cmdArgs[1]
    out_dir_mother <- cmdArgs[2]
    out_dir_father <- cmdArgs[3]
    
    mutTable <- read.table(input, sep = "\t", header = T, stringsAsFactors = F)
    mutTable <- mutTable[!is.na(mutTable$Phase_combined),]
    mutTable$context <- "NNNNN"
    
    for(proband in unique(mutTable$Proband_nr)) {
        
        output <- mutTable[mutTable$Proband_nr == proband,]
        
        output <- output[,c("Chr", "Pos_hg38", "Ref", "Alt", "context", "Proband_ref_count", "Proband_alt_count", "Phase_combined")]
        output_mother <- output[output$Phase_combined == "mother", -ncol(output)]
        output_father <- output[output$Phase_combined == "father", -ncol(output)]
        
        write.table(output_mother, file.path(out_dir_mother, paste0(proband, ".txt")), sep = "\t", quote = F, row.names = F, col.names = F)
        write.table(output_father, file.path(out_dir_father, paste0(proband, ".txt")), sep = "\t", quote = F, row.names = F, col.names = F)
        
    }
    
}

main()

main <- function(cmdArgs = commandArgs(T)) {
    
    overlap_map_file <- cmdArgs[1]
    out_randomazied_file <- cmdArgs[2]
    
    info <- file.info(overlap_map_file)
    
    if(any(info$size == 0)) {
        file.create(out_randomazied_file)
    } else {
        overlap_map <- read.table(overlap_map_file, sep = "\t", stringsAsFactors = F, header = F)
        from <- unique(overlap_map[,3])
        for(i in from) 
            overlap_map[ overlap_map[,3] == i, 4] <- sample(overlap_map[ overlap_map[,3] == i, 4])
        
        write.table(overlap_map, out_randomazied_file, sep = "\t", quote = F, col.names = F, row.names = F)
    }
    
}

main()

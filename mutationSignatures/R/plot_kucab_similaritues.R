library("dplyr")
library("tidyr")
library("gplots")
main <- function(cmdArgs = commandArgs(T)) {
    
    out_prefix <- cmdArgs[1]
    kucab_signatures_file <- cmdArgs[2]
    similarity_files <- cmdArgs[-(1:2)]
    
    out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/kucab_signatures/1_all"
    kucab_signatures_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/Mutagen53_sub_signature.txt"
    similarity_files <- list.files("/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/kucab_signatures/", pattern = "*kucab_similarity.txt", full = T)
    
    # Read files
    simi <- read_similarity_files(similarity_files)
    kucab_signatures <- read.table(kucab_signatures_file, check.names = F, sep = "\t", header = T)
    kucab_signatures <- setNames(colnames(kucab_signatures), paste0("V", 1:ncol(kucab_signatures)))
    
    
    # Convert similarities to two-dimensional matrix
    simi$tissue_sign <- paste0(simi$tissue, "_", simi$signature)
    simi_mat <- simi %>%
        select(tissue_sign, kucab_signature, cosine) %>%
        spread(kucab_signature, cosine) 
    
    # Appending names
    rownames(simi_mat) <- simi_mat[,1]
    simi_mat <- simi_mat[,-1]
    colnames(simi_mat) <- kucab_signatures[colnames(simi_mat)]
    simi_mat <- as.matrix(simi_mat)
    
    heatmap_cols <- colorRampPalette(c("#cee0ff", "#428aff"))(20)
    pdf(paste0(out_prefix, "_kucab_heatmap_cosine.pdf"))
    heatmap.2(simi_mat, Rowv = F, Colv = F, trace = "none", col = heatmap_cols, 
              density.info = "none", key.xlab = "Kucab et. al. signatures", symbreaks = F,
              colsep = 1:ncol(simi_mat), rowsep = 1:nrow (simi_mat), na.color = "grey",
              sepwidth=c(0.001,0.001), sepcolor = "black")
    dev.off()
    
    
}

read_similarity_files <- function(x) {
    
    results <- list()
    for(i in x) {
        current <- read.table(i, sep = "\t", stringsAsFactors = F, header = T)
        current[,1] <- gsub("(\\w+)-.+", "\\1", current[,1])
        colnames(current)[1] <- "tissue"
        colnames(current)[3] <- "kucab_signature"
        results[[i]] <- current
    }
    
    results <- do.call(rbind, results)
    rownames(results) <- 1:nrow(results)
    return(results)
}
                            

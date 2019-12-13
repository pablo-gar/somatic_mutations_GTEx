# assess mutation profile FCs similarity to mutations signatures
#


main <- function(cmdArgs=commandArgs(T)) {
    
    MUTATION_KEY <- c(`00` = "T>A", `01` = "T>T", `02` = "T>G", `03` = "T>C",
                      `10` = "C>A", `11` = "C>T", `12` = "C>G", `13` = "C>C")
    
    fc_file <- cmdArgs[1]
    signature_file <- cmdArgs[2]
    gene <- cmdArgs[3]
    signature_no <- cmdArgs[4]
    out_file <- cmdArgs[5]


    #fc_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations_targeted_mutProfile_differences/Lung/n6_0.0_0.7/1_profiles.txt"
    #signature_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cancerSignaturesStandard.txt"
    #gene <- "MLH1"
    #signature_no <- "Signature_6"
    #output_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations_targeted_mutProfile_differences_and_signatures/n6_0.0_0.7-Lung.txt"
    
    # Read signature file
    signatureTable <- read.table(signature_file, sep="\t", header=F, stringsAsFactors=F, row.names=1)
    colnames(signatureTable) <-paste0("Signature_", 1:ncol(signatureTable))
    
    mutationCode <- do.call(rbind, strsplit(rownames(signatureTable), "\\."))
    mut <- MUTATION_KEY[mutationCode[,1]]
    context <- mutationCode[,2]
    rownames(signatureTable) <- paste0(mut, '.', context)

    # Calculate frequencies
    signatureTable <- as.data.frame(t(t(signatureTable) / colSums(signatureTable)))
    
    
    # Read profiles
    fc_table <- read.table(fc_file, sep="\t", header=T, stringsAsFactors=F)
    fc_table <- fc_table[fc_table$gene == gene,]
    rownames(fc_table) <- paste0(fc_table$mut, ".", fc_table$context)
    muts <- unique(fc_table$mut)
    fc_table <- fc_table[,-(1:3)]
    
    # Gather profiles of desired signature and desired gene in the same order
    signatureCosine <- matrix(rep(signatureTable[,signature_no], ncol(fc_table)), byrow=F, ncol=ncol(fc_table))
    rownames(signatureCosine) <- rownames(signatureTable)
    fc_table <- as.matrix(log2(fc_table))
    
    cosine_results <- cosineSimilarity(signatureCosine, fc_table)
    
    results <- data.frame(mut=c("all", muts),
                          cosine=cosine_results[1], 
                          pvalue_right=sum(cosine_results[-1] >= cosine_results[1])/(length(cosine_results)-1),
                          pvalue_left=sum(cosine_results[-1] <= cosine_results[1])/(length(cosine_results)-1),
                          n_perm=length(cosine_results)-1,
                          stringsAsFactors=F
                          )
    
    for(i in 2:nrow(results)) {
        
        mut <- results$mut[i]
        
        current_sign <- signatureCosine[grep(mut, rownames(signatureCosine)),]
        current_fc <- fc_table[grep(mut, rownames(fc_table)),]
        cosine_results <- cosineSimilarity(current_sign[rownames(current_fc),], current_fc)
        
        results[i, c("cosine", "pvalue_right", "pvalue_left")] <- c(cosine_results[1],
                                                                    sum(cosine_results[-1] >= cosine_results[1])/(length(cosine_results)-1),
                                                                    sum(cosine_results[-1] <= cosine_results[1])/(length(cosine_results)-1)
                                                                    )
        
    }
    
    
    
    write.table(results, out_file, sep="\t", quote=F, col.names=T, row.names=F)
    
}

cosineSimilarity <- function(a, b) {
    
    stopifnot(dim(a) == dim(b))
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b))
    
    colSums(a * b) / (sqrt(colSums(a^2)) * sqrt(colSums(b^2)))
}

spearman_cor_matrix <- function(a,b) { 
    stopifnot(dim(a) == dim(b))
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b))
    
    results <- rep(0, ncol(a))
    
    for(i in 1:length(results))
        results[i] <- cor(a[,i], b[,i], method="spearman")
    
    return(results)
    
}

main()

# Looks at mutation profile differences between samples that have low/high expression
# of certain genes

#library("MatrixEQTL")
library("dplyr")
library("tidyr")
library("ggplot2")
#library("metap")
#source("../../R/FDR.R")
#source("../../R/plots.R")
source("../../R/gtex.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
#source("../../R/misc.R", chdir = T)
#source("../../R/mutationGeneAnnotation.R", chdir = T)

main <- function(cmdArgs=commandArgs(T)) {
    
    barColors <- c("#3BC9F3", "#2B2E34", "#FC3218", "#CAD0CE", "#9CD169", "#F1CAC9")
    
    HEIGHT = 2.2
    WIDTH = 4 
    NCOL_PLOT = 6 
    NROW_PLOT = 1
    
    nPerm = 1000
    

    expression_file <- cmdArgs[1]
    cntxt_file <- cmdArgs[2]
    output_prefix <- cmdArgs[3]


    #expression_file <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations_targeted_expressionMat/Lung/n6_0.0_0.7/expression_mat.txt"
    #cntxt_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationMatrix/Lung-n6_0.0_0.7.txt"
    #output_prefix <- "/scratch/users/paedugar/somaticMutationsProject/generalMutationAnalyses/results/mutation_expression_associations_targeted_mutProfile_differences/Lung/n6_0.0_0.7/"
    
    # Read files
    expression_mat <- read.table(expression_file, sep="\t", stringsAsFactors=F, header=T, check.names=F)
    cntxt_mat <- readMutations(cntxt_file, contextLength=3, ignoreN=T)
    colnames(cntxt_mat)[-(1:2)]<- sraToGtex(colnames(cntxt_mat)[-(1:2)], formatOut="long")
    
    # Calculate frequencies
    cntxt_mat[, -(1:2)] <- apply(cntxt_mat[, -(1:2)], 2, function(x) x/sum(x))
    
    # Calculate profile fc between individuals that have high/low expression for all genes
    quant <- 0.2
    results <- list()
    for(i in 1:nrow(expression_mat)) {
        
        # Getting top and bottom samples based on expression of given gene
        gene <- expression_mat[i, ncol(expression_mat)]
        exp_gene <- setNames(expression_mat[i, -ncol(expression_mat)], colnames(expression_mat)[-ncol(expression_mat)])
        exp_gene <- rankitNormalize(exp_gene)
        expression_quant <- quantile(exp_gene, probs = c(quant, 1-quant))
        
        top_samples <- names(exp_gene[exp_gene <= expression_quant[1]])
        bottom_samples <- names(exp_gene[exp_gene >= expression_quant[2]])
        
        # Create fc profiles of the averages across groups
        toPlot <- cntxt_mat[,1:2]
        toPlot$gene <- gene
        
        # Getting fc top/bottom
        top <- rowMeans(cntxt_mat[,colnames(cntxt_mat) %in% top_samples])
        bottom <- rowMeans(cntxt_mat[,colnames(cntxt_mat) %in% bottom_samples])
        toPlot$fc <- bottom / top
        
        # Getting fc of permutations
        index <- 3:ncol(cntxt_mat)
        top_l <- length(top)
        
        perms <- lapply(1:nPerm, function(x) {
                            topInd <- sample(index, top_l, replace=F)
                            bottomInd <- sample(index[!index %in% topInd], top_l, replace=T)
                            rowMeans(cntxt_mat[,topInd]) / rowMeans(cntxt_mat[,bottomInd])
                        })
        perms <- do.call(cbind, perms)
        colnames(perms) <- paste0('perm', 1:ncol(perms))
        
        # Putting together
        toPlot <- cbind(toPlot, perms)
        
        p <- ggplot(toPlot, aes(x=context, y=log2(fc), fill=mut)) +
                 geom_bar(stat='identity', position='dodge') +
                 facet_wrap(~mut, scale="free", nrow=NROW_PLOT) +
                 scale_fill_manual(values=barColors) + 
                 theme_grid_y() +
                 ylab(paste0("log2(fold-change profiles)\n(lowly/highly ", gene, " -expressing samples")) +
                 theme(axis.text.x=element_text(angle=45, hjust=1))
        
        ggsave(paste0(output_prefix, "gene_", gene, "_bottom_over_top_expressed_samples_profiles.pdf"), width = NCOL_PLOT * WIDTH, height = NROW_PLOT * HEIGHT)
        
        # Compile results
        results[[i]] <- toPlot
    }
    
    results <- do.call(rbind, results)
    write.table(results, paste0(output_prefix, "1_profiles.txt"), sep='\t', row.names=F, quote=F)
    
    
}

readMutations <- function(x, contextLength, ignoreN) {
        
    # Reads a file with columns as Mutation Type, Context, Ind 1, ... , Ind n
    # Merges rows with identical contexts depending on the context length
    
    MUTATION_KEY <- c(`00` = "T>A", `01` = "T>T", `02` = "T>G", `03` = "T>C",
                      `10` = "C>A", `11` = "C>T", `12` = "C>G", `13` = "C>C")

    #Getting number of columns
    con <- file(x, open = "r")
    columnsFile <- length(unlist(strsplit(readLines(con, n = 1), "\t")))
    close(con)


    x <- read.delim(x, header = T, sep = "\t", stringsAsFactors = F, colClasses = c("character", "character", rep("numeric", columnsFile - 2)) )

    if (ignoreN)
        x <- x[!grepl("N", x[,2]), ]

    currentLength <-  nchar(x[1,2])

    if (contextLength > currentLength)
        stop("context length specified greater than currently availabe")

    if (contextLength < currentLength) {
        lengthDiff <- (currentLength - contextLength) / 2 
        x[,2] <- substr(x[,2], lengthDiff + 1, currentLength - lengthDiff)

        # Put same context together
        x <- do.call(rbind, by(x, list(x[,1], x[,2]), function(x) {
                                   x[1,3:ncol(x)] <- colSums(x[,3:ncol(x)])
                                   return (x[1,])
                                    }   
                                )   
                    )   

    }   


    # Reordering and renaming
    x <- x[order(x[,1], x[,2]),]
    rownames(x) <- 1:nrow(x)
    x[,1] <- MUTATION_KEY[x[,1]]
    return(x)
}

rankitNormalize <- function(x) {
    qnorm((rank(x) - 0.5) / length(x))
}


main()

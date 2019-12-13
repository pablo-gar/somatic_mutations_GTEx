# Creates a scatter plot of shannon diversity vs number of mutations per sample
# Usage
# Rscript mutationVsTrasncriptomeDiversity.R mapFile1.txt [mapFile2.txt] [...]

source("../../R/gtex.R", chdir = T)
source("../../R/plots.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
#library("dplyr")

cmdArgs <- commandArgs(T)
outplot_TPM <- cmdArgs[1]
outplot_Reads <- cmdArgs[2]
mutationCountFiles <- cmdArgs[-(1:2)]

#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/Liver/n6_0.0_0.5/"
#mutationCountFiles <- file.path(workingDir, list.files(workingDir))

#----------------------------------#
# METHODS
#----------------------------------#
getMutationsPerGene <- function(x) {
    
    result <- by(x, x$sample, function(x) {
                 y <-as.data.frame(table(x$gene_id), stringsAsFactors = F)
                 colnames(y) <- c("gene_id", "mutations")
                 y$sample <- x$sample[1]
                 return(y)
                })
    return(do.call(rbind,result))
}


# calculates the shanon entropy of a vector
shannonEntropy <- function(x, percentage = T) {
    
    x <- x[ x != 0 & !is.na(x) ]
    
    shan <- -sum( (x/sum(x)) * log2(x/sum(x)) )
    
    if(percentage)
        shan <- shan / log2(length(x))
    
    return(shan)
}

# Calcautes the shanon entropy across columns in a matrix in an efficent way
# IT DOES NOT ELIMINATE 0s
shannonEntropyMatrix <- function(x) {
    
    small <- 0.000001
    
    # Doing frequencies by column
    x <- x / matrix(colSums(x), nrow = nrow(x), ncol = ncol(x), byrow=T)
    
    shan <- -colSums((x + small) * log2(x + small))
    
    return(shan)
}

appendShannon <- function(x, countFile, gtexId, colName, ignoreZero = T) {

    geneExp <- readAllGtexExpression(gtexId, countFile)
    
    if(ncol(geneExp) < 2 | nrow(geneExp) < 2)
        return(NULL)
    
    geneExp <- geneExp[,-1]
    
    if(ignoreZero) {
        shan <- apply(as.matrix(geneExp), 2, shannonEntropy)
    } else {
        shan <- shannonEntropyMatrix(as.matrix(geneExp))
    }
    

    x$new <- shan[x$gtexId]
    x <- x[!is.na(x$new),]
    colnames(x)[ncol(x)] <- colName
    
    return(x)
    
}

#' Randomizes the values in the columns of a data.frame
randomizeRows <- function(x, column) {
    x[,column] <- x[sample(1:nrow(x), nrow(x), replace = F),column]
    return(x)
}

#' Randomize the values in a column grouping by the "sample" column 
randomizePerSample <- function(mutationsPerGene, column) {
    
    
    results <- 
        mutationsPerGene %>%
        group_by(sample) %>%
        randomizeRows(column = column) %>%
        ungroup()
    
    return(results)
}

#----------------------------------#
# MAIN
#----------------------------------#

####
# Read mutations and process mutations
countTable <- data.frame(sample = "", mutations = rep(0, length(args)), stringsAsFactors = F)

for (i in 1:length(mutationCountFiles)) {
    mutation <- read.table(mutationCountFiles[i], sep = "\t", stringsAsFactors = F, header = F)
    countTable[i, "sample"] <- gsub(".txt", "", basename(mutationCountFiles[i]))
    countTable[i, "mutations"] <- sum(mutation)
}

# Adding gtexIds
gtexId <- setNames(sraToGtex(unique(countTable$sample), formatOut = "long"), unique(countTable$sample))
countTable$gtexId <- gtexId[countTable$sample]

#######
# Adding gene counts and TPM
countTable <- appendShannon(countTable,  CONFIG$auxiliaryFiles$expresionAllTissues, countTable$gtexId, "shannonDiversity_TPM")
countTable <- appendShannon(countTable,  CONFIG$auxiliaryFiles$readsGenesAllTissues, countTable$gtexId, "shannonDiversity_readCounts")

if(is.null(countTable)) {
    
    ggsave(outplot_TPM, ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100))
    ggsave(outplot_Reads, ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100))
    
} else {
    
    p1 <- scatter(countTable, x = "shannonDiversity_TPM", y = "mutations", regression = T)
    p2 <- scatter(countTable, x = "shannonDiversity_readCounts", y = "mutations", regression = T)

    ggsave(outplot_TPM, p1)
    ggsave(outplot_Reads, p2)
}

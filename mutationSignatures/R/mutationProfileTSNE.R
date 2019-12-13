#' Performs a tSNE analyses using only mutation profiles
#'
#' Author: Pablo E. Garcia-Nieto
#' Date: 3/4/2018
#'
#' Usage:
#'  Rscript mutatioProfileTSNE.R tSNE.pdf W1.txt [W2.txt] [...]
#'
#' @param integer context length of mutations, k-mer file has to be corresponding to this
#' @param filepath prefix for the output tSNE plot
#' @param filepath for a file with k-mer counts per gene, output of oligonucleotideGTF.R
#' @param filepaths [...] to W files to use for composite W
#'
#' @return Returns tSNE plots for gene expression, and mutation profile 

if(!dir.exists("../../R"))
       stop("Can't find the share R path, make sure to run script from directory where it lives")

library("Rtsne")
library("ggplot2")
library("reshape")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
source("./mutationSignatures/source/mutationSignaturesNMF.methods.R", chdir = T)

args <- commandArgs(T)
contextLength <- as.numeric(parseArg(args[1], sep = ","))
tsneOutputPrefix <- args[2]
oligoCountFile <- args[3]
Wfiles <- args[4:length(args)]


#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationMatrix/"
#Wfiles <- list.files(workingDir)[grep(".*n6.*.txt", list.files(workingDir))]
#Wfiles <- file.path(workingDir, Wfiles)
#contextLength <- 5
#tsneOutputPrefix <- "~/n6"
#oligoCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/GTF_gtex_oligoCounts_5.txt"

tsneOutputPrefix <- paste0(tsneOutputPrefix, "_context_", contextLength)
normilizeByExpression <- T
doPCA <- F
VAR_EXPLAINED_PCA <- 0.98


    


#--------------------------------------------
# METHODS

readWtables <- function(x, contextLength, ignoreN) {
    
    for(i in x) {
        W <- readMutations(i, contextLength = contextLength, ignoreN = T)
        rownames(W) <- paste0(W[,1], ".", W[,2])
        W <- W[, -(1:2)]
        colnames(W) <- paste0(basename(i), "|", colnames(W))
        if(exists("allWs")) {
            allWs <- cbind(allWs, W)
        } else {
            allWs <- W
        }
            
    }
    
    return(allWs)
}

#' Scatter plot for a tsne Obj with 2 dimensions
#' @param tsneObj - tsne matrix from the output from Rtsne package (x$Y)
#' @param categories - character, categories associated with each point in scatter
#' @param col - character, named vector where the names are the unique categories and the values their associated colors
plotTsne <- function(tsneObj, categories = NULL, col = tissueCol){
    
    require("ggplot2")
    stopifnot(ncol(tsneObj) == 2)
    
    return (plotColoredScatterFromMatrix(tsneObj, categories = categories, col = col))
}

#' Scatter plot from a 2-column matrix with colored points based on categories
#' @param x - matrix with 2 columns, the x and y coordinates
#' @param categories - character, categories associated with each point in scatter
#' @param col - character, named vector where the names are the unique categories and the values their associated colors

plotColoredScatterFromMatrix <- function(x, categories = NULL, col = tissueCol) {
    require("ggplot2")
    
    stopifnot(ncol(x) == 2)
    stopifnot(length(categories) == nrow(x))
    stopifnot(unique(sort(categories)) == unique(sort(names(col))))
    
    x <- as.data.frame(x, stringsAsFactors = F)
    
    x$categories <- categories
    x$col <- col[categories]
    
    shapes <- 1:length(col)
    names(shapes) <- names(col)
    x$shape <- shapes[categories]
    
    nPoints <- nrow(x)
    
    #x <- melt.data.frame(x, id.vars = c("categories", "col"), measure.vars = 1:2)
    
    p <- ggplot(x, aes(x = V1, y = V2)) +
    theme_noGrid()
    
    if(is.null(categories)) {
        p <- p + geom_point()
    } else {
        
        p <- p + geom_point( colour = x$col, shape = x$shape, size = 2.5)
        
        
        forLegend <- data.frame(x = 1, y = length(col):1, label = names(col), col = col, shape = shapes[names(col)])
        pLegend <-  ggplot(forLegend, aes(x = x, y = y)) +
            geom_point(colour = forLegend$col, shape = forLegend$shape, size = 6) +
            geom_label(aes(x = x + 0.05, label = label), hjust = 0 ,colour = "white", fill = forLegend$col, fontface = "bold", size = 3) +
            xlim(0.5, 2.5) +
            ylim(-2, nrow(forLegend) + 3) +
            theme_void() 
    }
    
    
    return(list(scatter = p, legend = pLegend))
}

performTSNE <- function(x, doPCA = T, minVarExplained = 0.9) {
    
    if(doPCA) { 
        pcaResult <- prcomp(t(x))
        varExplained <- cumsum(pcaResult$sdev^2 / sum(pcaResult$sdev^2))
        if(all(!(varExplained <= minVarExplained))) {
            x <- pcaResult$rotation
        } else {
            x <- pcaResult$rotation[, varExplained <= minVarExplained]
        }
    }
    
    #Eliminate duplicated columns
    x <- x[, !duplicated(colnames(x))]
    
    tsne <- Rtsne(x, dims = 2, max_iter = 500, verbose = T)
    
    return(tsne)
    
}

performPCA <- function(x) {
    
    #Eliminate duplicated columns
    x <- x[, !duplicated(colnames(x))]
    pcaResult <- prcomp(t(x))
    
    return(pcaResult)
    
}
    

readGeneExpression <- function(sraIds) { 
    
    gtexIds <- sraToGtex(sraIds, formatOut = "long")
    geneExp <- readAllGtexExpression(gtexIds)
    rownames(geneExp) <- geneExp[,1]
    geneExp <- geneExp[,-1]
    colnames(geneExp) <- gtexLongToSra(colnames(geneExp))
    
    return(geneExp)
}

#' Reads the kmer count file x
readOligoCountPerGene <- function(x) {
    x <- read.table(x, sep = "\t", header = T, stringsAsFactors = F)
    rownames(x) <- x[,1]
    x <- x[,-1]
}


#' Normilizes mutation count by sum(contextCount*expressionOfGene), sum is over all genes
#' @param mutationTable data.frame - with columns TODO
#' @param oligoCountFile path to file with kmer counts, columns gene_id, kmer_1, kmer_2, ..., kmer. K should be the same as contextLenght of mutation
normalizeByExpression <- function(mutationTable, oligoCountFile, geneExp) {
    
    
    oligoCount <- readOligoCountPerGene(oligoCountFile)
    
    context <- gsub(".+\\.(.+)", "\\1", rownames(mutationTable))
    
    if(nchar(context[1]) != nchar(colnames(oligoCount)[1]))
        stop("Context are not compatible between the muation Table and the k-mer table")
    
    # Loading gene Expression
    sraIds <- gsub(".+\\|(.+)", "\\1", colnames(mutationTable))
    geneExp <- geneExp[rownames(oligoCount), ]
    
    # Calculating total oligocount in transcriptome of each individual 
    oligoCountTotal <-  as.matrix(t(oligoCount)) %*% as.matrix(geneExp)
    
    # Selecting the same individuals
    #colnames(oligoCountTotal) <- gtexLongToSra(colnames(oligoCountTotal))
    mutationTable <- mutationTable[ , sraIds %in% colnames(oligoCountTotal)]
    sraIds <- gsub(".+\\|(.+)", "\\1", colnames(mutationTable))
    oligoCountTotal <- oligoCountTotal[, sraIds]
    colnames(oligoCountTotal) <- colnames(mutationTable)
    
    # Compacting context counts, e.g. ACA = ACA + AGA
    middleBaseI <- ceiling(nchar(rownames(oligoCountTotal)[1]) / 2)
    bases <- c("T", "T", "C", "C")
    names(bases) <- c("A", "T", "G", "C")
    compactContext <- paste0(substring(rownames(oligoCountTotal), 1, middleBaseI - 1), 
                             bases[substring(rownames(oligoCountTotal), middleBaseI, middleBaseI)] , 
                             substring(rownames(oligoCountTotal), middleBaseI + 1,))
    compactOligoCount <- matrix(0, ncol = ncol(oligoCountTotal), nrow = nrow(oligoCountTotal) / 2 , dimnames = list(unique(compactContext), colnames(oligoCountTotal)))
    for(i in 1:nrow(compactOligoCount)) { 
        currentContext <- rownames(compactOligoCount)[i]
        compactOligoCount[i,] <- colSums(oligoCountTotal[compactContext == currentContext, ])
    }
    
    # Returning the table
    compactOligoCount <- compactOligoCount[context,]
    
    return(compactOligoCount)
}

doSaveTSNE <- function(Wjoint, doPCA, minVarExplained, tsneOutputPrefix, useFreq = T) {
    
    # Calculate freqs
    if(useFreq)
        Wjoint <- as.data.frame(t(t(Wjoint) / colSums(Wjoint)))
    
    Wjoint <- t(as.matrix(Wjoint))
    Wjoint <- Wjoint[!duplicated(Wjoint),]
    
    
    #############
    # Do tsne
    tsne <- performTSNE(Wjoint, doPCA = doPCA, minVarExplained = minVarExplained)
    
    # Get colors for tissues
    tissues <- gsub("(.+)\\|.+", "\\1", rownames(Wjoint))
    tissueCol <- rainbow(length(unique(tissues)))
    names(tissueCol) <- unique(sort(tissues))
    
    # Plot tsne
    p <- plotTsne(tsne$Y, categories = tissues, col = tissueCol)
    
    ggsave(paste0(tsneOutputPrefix, "scatterMut.pdf"), p[[1]])
    ggsave(paste0(tsneOutputPrefix, "legend.pdf"), p[[2]])
    
    #############
    # Do PCA
    pcaResult <- performPCA(Wjoint)
    pcaVarExplained <- (pcaResult$sdev^2 / sum(pcaResult$sdev^2)) * 100
    pcaResult <- pcaResult$rotation[, 1:2]
    
    # For some reason it is good to eliminate highest and lowest value on PC2
    pcaResult <- pcaResult[ pcaResult[,2] < max(pcaResult[,2]) & pcaResult[,2] > min(pcaResult[,2]), ]
    colnames(pcaResult) <- c("V1", "V2")
    tissues <- gsub("(.+)\\|.+", "\\1", rownames(pcaResult))
    tissueCol <- rainbow(length(unique(tissues)))
    names(tissueCol) <- unique(sort(tissues))
    
    # Plot PCA
    p2 <- plotTsne(pcaResult, categories = tissues, col = tissueCol)
    p2[[1]] <- p2[[1]] + xlab(paste0("PC1 (", round(pcaVarExplained[1]), "% var explained)")) + ylab(paste0("PC2 (", round(pcaVarExplained[2]), "% var explained)")) 
    
    ggsave(paste0(tsneOutputPrefix, "PCA_scatterMut.pdf"), p2[[1]])
    ggsave(paste0(tsneOutputPrefix, "PCA_legend.pdf"), p2[[2]])
    
}

#--------------------------------------------

#--------------------------------------------
# MAIN


#-------
# Work with mutations
Wjoint <- readWtables(Wfiles, contextLength)
# Exclude exo samples
Wjoint <- Wjoint[, !grepl("EXO", colnames(Wjoint))]

sraIds <- gsub(".+\\|(.+)", "\\1", colnames(Wjoint))



# Normalize by freqeuncy of nucleotides in entire transcriptome
if(normilizeByExpression) {
    geneExp <- readGeneExpression(sraIds)
    contextCounts <- normalizeByExpression(Wjoint, oligoCountFile, geneExp)
}

# Perform and save tSNE
doSaveTSNE(Wjoint = Wjoint, doPCA = doPCA, minVarExplained = VAR_EXPLAINED_PCA, tsneOutputPrefix = paste0(tsneOutputPrefix, "_original_"))

if(normilizeByExpression) {
    Wjoint_copy <- Wjoint
    Wjoint_copy[,colnames(contextCounts)] <- Wjoint_copy[,colnames(contextCounts)] / contextCounts
    
    doSaveTSNE(Wjoint = Wjoint_copy, doPCA = doPCA, minVarExplained = VAR_EXPLAINED_PCA, tsneOutputPrefix = paste0(tsneOutputPrefix, "_normalized_"))
    doSaveTSNE(Wjoint = contextCounts, doPCA = doPCA, minVarExplained = VAR_EXPLAINED_PCA, tsneOutputPrefix = paste0(tsneOutputPrefix, "_countsOnly_"))
}

# DONE working with mutations
#-----------

##------
## Work with expression
#
## Reading expresion
sraIds <- gsub(".+\\|(.+)", "\\1", colnames(Wjoint))

if(!normilizeByExpression)
    geneExp <- readGeneExpression(sraIds)

Wjoint <- Wjoint[ , sraIds %in% colnames(geneExp)]
sraIds <- gsub(".+\\|(.+)", "\\1", colnames(Wjoint))
geneExp <- geneExp[,sraIds]
colnames(geneExp) <- colnames(Wjoint)

doSaveTSNE(Wjoint = geneExp, doPCA = doPCA, minVarExplained = VAR_EXPLAINED_PCA, tsneOutputPrefix = paste0(tsneOutputPrefix, "_geneExpression_"), useFreq = F)
#
## Getting tissues ready
#tissues <- gtexLongToTissues(rownames(geneExp))
#tissueCol = rainbow(length(unique(tissues)))
#names(tissueCol) = unique(sort(tissues))
#
## Plot tsne
#p <- plotTsne(tsneExpression, categories = tissues, col = tissueCol)
#
#prefix <- paste0(tsneOutputPrefix, "tsne-expression_")
#ggsave(paste0(prefix, "scatterMut.pdf"), p[[1]])
#ggsave(paste0(prefix, "legend.pdf"), p[[2]])
#
##--------------------------------------------

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
library("cluster")
source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
source("../../R/plots.R", chdir = T)
source("../../mutationSignatures/R/mutationSignatures/source/mutationSignaturesNMF.methods.R", chdir = T)

args <- commandArgs(T)
contextLength <- as.numeric(parseArg(args[1], sep = ","))
tsneOutputPrefix <- args[2]
oligoCountFile <- args[3]
readCountFile <- args[4]
transDiversityFile <- args[5]
Wfiles <- args[6:length(args)]

color_point <- rbind(
                     data.frame(a = "Adipose_Subcutaneous", b = "#e6194B"),
                     data.frame(a = "Adipose_Visceral_Omentum", b = "#b51255"),
                     data.frame(a = "Adrenal_Gland", b = "#ffe119"),
                     data.frame(a = "Artery_Aorta", b = "#f58231"),
                     data.frame(a = "Artery_Coronary", b = "#4363d8"),
                     data.frame(a = "Artery_Tibial", b = "#911eb4"),
                     data.frame(a = "Brain_Caudate_basal_ganglia", b = "#1d4d00"),
                     data.frame(a = "Brain_Cerebellar_Hemisphere", b = "#276600"),
                     data.frame(a = "Brain_Cerebellum", b = "#318000"),
                     data.frame(a = "Brain_Cortex", b = "#3b9900"),
                     data.frame(a = "Brain_Frontal_Cortex_BA9", b = "#44b300"),
                     data.frame(a = "Brain_Hippocampus", b = "#50d100"),
                     data.frame(a = "Brain_Hypothalamus", b = "#58e600"),
                     data.frame(a = "Brain_Nucleus_accumbens_basal_ganglia", b = "#62ff00"),
                     data.frame(a = "Brain_Putamen_basal_ganglia", b = "#71ff1a"),
                     data.frame(a = "Breast_Mammary_Tissue", b = "#7c7c7c"),
                     data.frame(a = "Colon_Sigmoid", b = "#42d4f4"),
                     data.frame(a = "Colon_Transverse", b = "#f032e6"),
                     data.frame(a = "Esophagus_Gastroesophageal_Junction", b = "#bfef45"),
                     data.frame(a = "Esophagus_Mucosa", b = "#469990"),
                     data.frame(a = "Esophagus_Muscularis", b = "#fabebe"),
                     data.frame(a = "Heart_Atrial_Appendage", b = "#e28d2b"),
                     data.frame(a = "Heart_Left_Ventricle", b = "#c67d29"),
                     data.frame(a = "Liver", b = "#800000"),
                     data.frame(a = "Lung", b = "#aaffc3"),
                     data.frame(a = "Muscle_Skeletal", b = "#808000"),
                     data.frame(a = "Nerve_Tibial", b = "#ffd8b1"),
                     data.frame(a = "Ovary", b = "#000075"),
                     data.frame(a = "Pancreas", b = "#a9a9a9"),
                     data.frame(a = "Pituitary", b = "#000000"),
                     data.frame(a = "Prostate", b = "#9faa00"),
                     data.frame(a = "Skin_Not_Sun_Exposed_Suprapubic", b = "#ff9868"),
                     data.frame(a = "Skin_Sun_Exposed_Lower_leg", b = "#ff8147"),
                     data.frame(a = "Small_Intestine_Terminal_Ileum", b = "#681111"),
                     data.frame(a = "Spleen", b = "#757f68"),
                     data.frame(a = "Stomach", b = "#9b78bc"),
                     data.frame(a = "Testis", b = "#7c7c7c"),
                     data.frame(a = "Thyroid", b = "#9db4dd"),
                     data.frame(a = "Whole_Blood", b = "#ff0000"),
                     data.frame(a = "Whole_Blood_EXO", b = "#2b2b2b")
                     )
color_point <- apply(color_point, 2, as.character)

#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationMatrix/"
#Wfiles <- list.files(workingDir)[grep(".*n6.*.txt", list.files(workingDir))]
#Wfiles <- file.path(workingDir, Wfiles)
#Wfiles <- Wfiles[grep("n6_0.0_0.7", Wfiles)]
#Wfiles <- Wfiles[!grepl("EXO", Wfiles)]
#contextLength <- 5
#tsneOutputPrefix <- "~/n6"
#oligoCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/GTF_gtex_oligoCounts_5.txt"
#readCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_OLD"
#transDiversityFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/transcriptome_shannon_diversity.txt"

tsneOutputPrefix <- paste0(tsneOutputPrefix, "_context_", contextLength)
normilizeByExpression <- T
doPCA <- F
VAR_EXPLAINED_PCA <- 0.98


    


#--------------------------------------------
# METHODS

readWtables <- function(x, contextLength, ignoreN) {
    for(i in x) {
        flush.console()
        cat(i, "\n")
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
plotTsne <- function(tsneObj, categories = NULL, col = tissueCol, extraCol = NULL, alpha = NULL){
    
    require("ggplot2")
    stopifnot(ncol(tsneObj) == 2)
    
    return (plotColoredScatterFromMatrix(tsneObj, categories = categories, col = col, extraCol = extraCol, alpha = alpha))
}

#' Scatter plot from a 2-column matrix with colored points based on categories
#' @param x - matrix with 2 columns, the x and y coordinates
#' @param categories - character, categories associated with each point in scatter
#' @param col - character, named vector where the names are the unique categories and the values their associated colors

plotColoredScatterFromMatrix <- function(x, categories = NULL, col = tissueCol, extraCol = NULL, alpha = NULL) {
    require("ggplot2")
    
    stopifnot(ncol(x) == 2)
    stopifnot(length(categories) == nrow(x))
    stopifnot(unique(sort(categories)) == unique(sort(names(col))))
    
    x <- as.data.frame(x, stringsAsFactors = F)
    x$categories <- categories
    
    # Creating col column
    if(is.null(extraCol)) {
        x$col <- col[categories]
    } else {
        x$col <- extraCol
    }
    
    # Creating alpha column
    if(!is.null(alpha))
        x$alpha <- alpha
    
    shapes <- rep(1:19, length.out = length(col))
    names(shapes) <- names(col)
    x$shape <- shapes[categories]
    
    nPoints <- nrow(x)
    
    #x <- melt.data.frame(x, id.vars = c("categories", "col"), measure.vars = 1:2)
    p <- ggplot(x, aes(x = V1, y = V2)) +
        xlab("tSNE dimension 1") +
        ylab("tSNE dimension 2") +
        theme_noGrid()
    
    if(is.null(categories)) {
        p <- p + geom_point()
    } else {
        
        if (is.null(extraCol)) {
            if(is.null(alpha)) {
                p <- p + geom_point(colour = x$col, shape = x$shape, size = 1.5)
            } else {
                p <- p + geom_point(alpha = alpha, colour = x$col, shape = x$shape, size = 1.5)
            }
        } else {
            p <- p + geom_point(aes(colour = col), shape = x$shape, size = 1.5) + scale_colour_distiller(palette = "Spectral")
        }
            
        #p <- p + geom_point(colour = x$col, size = 1.0)
        
        
        forLegend <- data.frame(x = 1, y = length(col):1, label = names(col), col = col, shape = shapes[names(col)])
        pLegend <-  ggplot(forLegend, aes(x = x, y = y)) +
            geom_point(colour = forLegend$col, shape = forLegend$shape, size = 6) +
            #geom_point(colour = forLegend$col, size = 6) +
            geom_label(aes(x = x + 0.05, label = label), hjust = 0 ,colour = "white", fill = forLegend$col, fontface = "bold", size = 3) +
            xlim(0.5, 2.5) +
            ylim(-2, nrow(forLegend) + 3) +
            theme_void() 
    }
    
    
    return(list(scatter = p, legend = pLegend))
}

performTSNE <- function(x, doPCA = T, minVarExplained = 0.9, seed = 1) {
    
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
    
    set.seed(seed)
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

doSaveTSNE <- function(Wjoint, doPCA, minVarExplained, tsneOutputPrefix, useFreq = T, do_extended_plots = F, alphas = NULL, tissues_sil = NULL) {
    
    # Calculate freqs
    if(useFreq)
        Wjoint <- as.data.frame(t(t(Wjoint) / colSums(Wjoint)))
    
    Wjoint <- t(as.matrix(Wjoint))
    Wjoint <- Wjoint[!duplicated(Wjoint),]
    
    #############
    # Do tsne
    tsne <- performTSNE(Wjoint, doPCA = doPCA, minVarExplained = minVarExplained, seed = 2)
    
    # Get colors for tissues
    tissues <- gsub("(.+)-.+\\|.+", "\\1", rownames(Wjoint))
    sraIds <- gsub(".+-.+\\|(.+)", "\\1", rownames(Wjoint))
    names(tissues) <- sraIds
          
    categories <- color_point[ color_point[,1] %in% tissues, 1]
    tissueCol <- color_point[ color_point[,1] %in% tissues, 2]
    names(tissueCol) <- categories
    
    # Plot tsne
    p <- plotTsne(tsne$Y, categories = tissues, col = tissueCol)
    ggsave(paste0(tsneOutputPrefix, "scatterMut.pdf"), p[[1]])
    ggsave(paste0(tsneOutputPrefix, "legend.pdf"), p[[2]], height = 12)
    
    # DO more plots: focusing on tissues, metadata, sihloutte
    if(do_extended_plots) {
        
        #########
        # WORKING with alphas i.e. higlighting some tissues
        if(!is.null(alphas)) {
            for(tissue_focus in alphas) {
                alpha <- rep(0.07, length(tissues))
                alpha[grep(tissue_focus, tissues)] <- 1
                p <- plotTsne(tsne$Y, categories = tissues, col = tissueCol, alpha = alpha)
                ggsave(paste0(tsneOutputPrefix, "scatterMut_alpha_", tissue_focus, ".pdf"), p[[1]])
            }
        }
            
        #DONE
        #########
        
        ########
        # Doing shiloutte
        dist_mat <- dist(tsne$Y)
        
        # Get individual widths
        results <- list()
        random_avg <- vector(mode = "numeric")
        # groups to try, user-defined clusters, all tissues and 10 random people
        nPeople <- 20
        tissues_sil <- c(tissues_sil, unique(tissues), sraIds[sample(length(sraIds), nPeople, replace = F)])
        for(i in seq_along(tissues_sil)) {
            
            current_tissue <- tissues_sil[i]
            
            # Real clusters
            if(grepl("SRR", current_tissue)) {
                clusters <- setNames(rep(1, length(unique(sraIds))), unique(sraIds))
                label <- sraIds
            } else {
                clusters <- setNames(rep(1, length(unique(tissues))), unique(tissues))
                label <- tissues
            }
            
            clusters[grep(current_tissue, names(clusters))] <- 2
            
            sil <- silhouette(clusters[label], dist_mat)
            avg_widths_conf_int <- bootstrap_confidence_interval(sil[sil[,1] == 2,3], mean, bootstrap_counts = 1e4)
            results[[i]] <- data.frame(avg = mean(sil[sil[,1] == 2 ,3]), conf_neg = avg_widths_conf_int[1], conf_pos = avg_widths_conf_int[2], tissue = current_tissue, stringsAsFactors = F)
            
            # Random clusters
            for(i in 1:10) {
                # Randomize tissue labels 
                sil <- silhouette(clusters[label[sample(length(label))]], dist_mat)
                random_avg <- c(random_avg, mean(sil[sil[,1] == 2 ,3]))
            }
            
        }
        
        # compile widths
        results <- as.data.frame(do.call(rbind, results))
        # compile people into single point
        people <- results[grepl("SRR", results$tissue),]
        results <- results[!grepl("SRR", results$tissue),]
        results <- rbind(results, data.frame(avg = mean(people$avg), conf_neg = mean(people$conf_neg), conf_pos = mean(people$conf_pos), tissue = "People average", stringsAsFactors = F))
        
        # plot
        results <- results[order(-results$avg),]
        results$tissue <- factor(results$tissue, levels = results$tissue, ordered = T)
        p_sil <- pointRange(results, x = "tissue", y = "avg", errorBarYmax = "conf_pos", errorBarYmin = "conf_neg", pSize = 0.5)
        p_sil <- p_sil + 
            geom_hline(yintercept = mean(random_avg), linetype = "dashed", colour = "lightskyblue3") +
            annotate("rect", fill = "lightskyblue3", alpha = 0.5, xmin = -Inf, xmax = Inf, ymin = mean(random_avg) - 2*sd(random_avg), ymax = mean(random_avg) + 2*sd(random_avg)) +
            theme_grid_y() +
            ylab("Avergae silhouette") +
            xlab("") +
            theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
        ggsave(paste0(tsneOutputPrefix, "scatterMut_silhouette.pdf"), p_sil, width = 10, height = 4)
        write.table(results, paste0(tsneOutputPrefix, "scatterMut_silhouette.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
        #DONE
        ## 
        
        
        #########
        # WORKING with METDATA
        # Read metadata
        rownames(tsne$Y) <- sraIds
        tsne_meta <- read_all_metadata(sraIds, readCountFile, transDiversityFile)
        shared_inds <- sraIds[sraIds %in% rownames(tsne_meta)]
        tsne_meta <- cbind(as.data.frame(tsne$Y[shared_inds,]), tsne_meta[shared_inds,])
        tissues <- tissues[shared_inds]
        categories <- color_point[ color_point[,1] %in% tissues, 1]
        tissueCol <- color_point[ color_point[,1] %in% tissues, 2]
        names(tissueCol) <- categories
        for(feature in colnames(tsne_meta)[-(1:2)]) {
            p <- plotTsne(tsne_meta[,1:2], categories = tissues, col = tissueCol, extraCol = tsne_meta[,feature])
            ggsave(paste0(tsneOutputPrefix, "scatterMut_metadata_", feature, ".pdf"), p[[1]])
        }
        # DONE
        ########
    }
    
    
    ##############
    ## Do PCA
    #pcaResult <- performPCA(Wjoint)
    #pcaVarExplained <- (pcaResult$sdev^2 / sum(pcaResult$sdev^2)) * 100
    #pcaResult <- pcaResult$rotation[, 1:2]
    #
    ## For some reason it is good to eliminate highest and lowest value on PC2
    ##pcaResult <- pcaResult[ pcaResult[,2] < max(pcaResult[,2]) & pcaResult[,2] > min(pcaResult[,2]), ]
    #
    #colnames(pcaResult) <- c("V1", "V2")
    #tissues <- gsub("(.+)\\|.+", "\\1", rownames(pcaResult))
    #tissueCol <- rainbow(length(unique(tissues)))
    #names(tissueCol) <- unique(sort(tissues))
    #
    ## Plot PCA
    #p2 <- plotTsne(pcaResult, categories = tissues, col = tissueCol)
    #p2[[1]] <- p2[[1]] + xlab(paste0("PC1 (", round(pcaVarExplained[1]), "% var explained)")) + ylab(paste0("PC2 (", round(pcaVarExplained[2]), "% var explained)")) 
    #
    #ggsave(paste0(tsneOutputPrefix, "PCA_scatterMut.pdf"), p2[[1]])
    #ggsave(paste0(tsneOutputPrefix, "PCA_legend.pdf"), p2[[2]])
    
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
doSaveTSNE(Wjoint = Wjoint, doPCA = doPCA, useFreq = F, minVarExplained = VAR_EXPLAINED_PCA, tsneOutputPrefix = paste0(tsneOutputPrefix, "_original_"))

if(normilizeByExpression) {
    Wjoint_copy <- Wjoint
    Wjoint_copy[,colnames(contextCounts)] <- Wjoint_copy[,colnames(contextCounts)] / contextCounts
    
    doSaveTSNE(Wjoint = Wjoint_copy, doPCA = doPCA, useFreq = F, minVarExplained = VAR_EXPLAINED_PCA, tsneOutputPrefix = paste0(tsneOutputPrefix, "_normalized_"), 
               do_extended_plots = T, alphas = c("Skin", "Spleen", "Lung", "Liver", "Esophagus_Muco", "Whole_Blood", "Thyroid", "Artery_Coronary", "Stoma|Colon_T|Small", "Muscle|Heart"),
               tissues_sil = c("Stoma|Colon_T|Small", "Brain", "Skin", "Muscle|Heart"))
    doSaveTSNE(Wjoint = contextCounts, doPCA = doPCA, useFreq = F, minVarExplained = VAR_EXPLAINED_PCA, tsneOutputPrefix = paste0(tsneOutputPrefix, "_countsOnly_"))
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

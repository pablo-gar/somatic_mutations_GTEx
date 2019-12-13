## Given a series of mutation signature table (tables W from NMF) 
## Performs hierarchical clustering and then k-means clustering for a define set of clusters
## The number of clusters shoud be picked manually based on the hierarchical clustering
##
## Author: Pablo E. Garcia-Nieto
## Date: 2/28/2018
##
## Usage:
##  Rscript singatureCancerSimilarity.R n6_ 0.01 W1.txt [W2.txt] [...]
##
## @arg filepath prefix for all output files, this include:
##      - [prefix]_heatmap.pdf
##      - [prefix]_dendrogram.pdf
##      - [prefix]_consensusSignatures.txt
##      - [prefix]_correspondanceSignatures.txt
## @arg filepath for output heatmap plot that has the clustering of mutation matrices
## @arg intger, heigh to which cut hte dendrogram
## @arg filepaths [...] to W files to use for composite W
##
## @return Returns composite W and heatmap with hierarchical clustering

library("gplots")
#library("wesanderson")

args <- commandArgs(T)
outPrefix <- args[1]
heightCut <- as.numeric(args[2])
Wfiles <- args[3:length(args)]

#outPrefix <- "n6"
#heightCut <- 0.15
#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/"
#Wfiles <- list.files(workingDir)[grep(".*n6.*W.txt", list.files(workingDir))]
#Wfiles <- file.path(workingDir, Wfiles)
               

MUTATION_COLORS <- c(`C>A`="#3BC9F3", `C>G` = "#2B2E34", `C>T` = "#FC3218", `T>A` = "#CAD0CE", `T>C` = "#9CD169", `T>G` = "#F1CAC9")
SIGNATURE_COLORS <- c("gray73" , "gray43")#wes_palette(n = 2, name="GrandBudapest")
HEATMAP_GRADIENT <- colorRampPalette(c("#FFF5F0", "#A50F15"))(20)
HEATMAP_CELL_SIZE <- 0.2

#--------------------------------------------
# METHODS

readWtables <- function(x) {
    
    for(i in x) {
        W <- read.table(i, sep = "\t", header = F, row.names = 1, check.names = F)
        colnames(W) <- paste0(basename(i), "|", colnames(W))
        if(exists("allWs")) {
            allWs <- cbind(allWs, W)
        } else {
            allWs <- W
        }
            
    }
    
    return(allWs)
}


createHeatmap <- function(x, cutHeight = NULL, transpose = T, heatmapFile) {
    require(gplots)
    x <- as.matrix(x)
    
    #------
    # Get all labels ready
    numericMut <- gsub("(\\d+)\\.\\w+", "\\1", rownames(x))
    contextMut <- gsub("\\d+\\.(\\w+)", "\\1", rownames(x))
    stringMut <- toStringMutation(numericMut)
    
    #Order matrix based on mutation types
    x <- x[order(stringMut, contextMut), ]
    stringMut <- sort(stringMut)
    
    #Get color bar for mutation types
    mutColorBar <- MUTATION_COLORS[stringMut]
    #-----
    
    
    #-----
    # H clustering
    hc <- hclust(dist(t(x)))
    hcDen <- as.dendrogram(hc)
    
    if(!is.null(cutHeight)){
        membership <- cutree(hc, h = cutHeight)[hc$order]
        # Getting signture at bottom as 1, then 2, ....
        tempNames <- names(membership)
        membership <- mapVectors(membership, 1:max(membership)) 
        names(membership) <- tempNames
    }
    #-----
    
    
    #----
    # Gathering pars for heatmap2
    heatPars <- list(x = x,
                     col = HEATMAP_GRADIENT, 
                     dendrogram = "col", Rowv = FALSE, Colv = hcDen,
                     symbreaks = F,
                     RowSideColors = mutColorBar,
                     #scale = "column",
                     trace = "none", density = "none"
                     )
    heightPlot <- nrow(x) * HEATMAP_CELL_SIZE
    widthPlot <- ncol(x) * HEATMAP_CELL_SIZE
    
    if(!is.null(cutHeight)){
        colors <- mapVectors(membership, SIGNATURE_COLORS)
        names(colors) <- names(membership)
        colors <- colors[colnames(x)]
        heatPars <- c(heatPars, list(ColSideColors = colors))
        
        if(transpose) {
            heatPars$RowSideColors = colors
        }
    }
    
    if(transpose) {
        heatPars$x <- t(x)
        heatPars$dendrogram = "row"
        heatPars$Rowv = hcDen
        heatPars$Colv = FALSE
        heatPars$ColSideColors <- mutColorBar
        
        heightPlot <- ncol(x) * HEATMAP_CELL_SIZE
        widthPlot <- nrow(x) * HEATMAP_CELL_SIZE
        
    }
    
    #----
    
    pdf(heatmapFile, height = heightPlot, width = widthPlot)
    p <- do.call(heatmap.2, heatPars)
    dev.off()
    
    return(list(heatmap = p, hc = hc, membership = membership))
    
}

#' Converts mutations from numeric format to stringFormat
toStringMutation <- function(x) {
    MUTATION_KEY <- c(`00` = "T>A", `01` = "T>T", `02` = "T>G", `03` = "T>C",
                      `10` = "C>A", `11` = "C>T", `12` = "C>G", `13` = "C>C")
    return(MUTATION_KEY[x])
}

#' Returns the best number of clusters based on minimizing within clusters sum of squares
getBestKClusters <- function(x) {
    
    require(NbClust)
    require(factoextra)
    
    p <- fviz_nbclust(x, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)
    
    return(p)
}

#' Gets a list of intercalating colors based on a previous series of color
expandColors <- function(x, n) {
    
    sizeX <- length(x)
    stopifnot(n > sizeX)
    
    result <- vector("character", n)
    lastIndex <- sizeX*(floor(n/sizeX))
    result[1:lastIndex] <- rep(x, floor(n/sizeX))
    if(n%%sizeX > 0)
        result[(lastIndex + 1): length(result)] <- x[1:(n%%sizeX)]
    
    return(result)
    
}

#' Maps any vector to a numeric vector, whenever there is a change
#' in x in maps the next element in y
mapVectors <- function(x, y) {
    
    result <- vector(mode(y), length(x))
    j <- 1
    result[1] <- y[j]
    
    previous  <- x[1]
    for(i in 2:length(x)) {
        
        if(previous != x[i]) 
            j <- ifelse(j+1 > length(y), 1, j+1)
        
        result [i] <- y[j]
        previous <- x[i]
    }
    
    return(result)
}
    
#' Averages colums belonging to same groups defined in groups
#' @arg x - data.frame of matrix
#' @arg groups - named vectors containing integer as the groups and the names correspond to colnames of x
colAverageByGroups <- function(x, groups) {
    
    result <- matrix(0, ncol = length(unique(groups)), nrow = nrow(x))
    rownames(result) <- rownames(x)
    
    uniqGroups <- groups[!duplicated(groups)]
    
    for(i in 1:length(uniqGroups)) {
        belonging <- names(groups)[groups == uniqGroups[i]]
        result[,i] <- rowMeans(as.matrix(x[, colnames(x) %in% belonging]))
    }
    
    return(result)
}
        
        

#--------------------------------------------

#--------------------------------------------
# MAIN


Wjoint <- readWtables(Wfiles)

# Calculate freqs
Wjoint <- as.data.frame(t(t(Wjoint) / colSums(Wjoint)))

# Get heatmap
heatmapResult <- createHeatmap(Wjoint, heightCut, transpose = T, heatmapFile = paste0(outPrefix, "_heatmap.pdf"))
pdf(paste0(outPrefix, "_dendrogram.pdf"))
plot(heatmapResult[[2]], hang = -1)
abline(h = heightCut)
dev.off()

# Create consensus signatures
consensusSignatures <- colAverageByGroups(Wjoint, groups = heatmapResult$membership)
write.table(consensusSignatures, paste0(outPrefix, "_consensusSignatures.txt"), col.names = F, row.names = T, sep = "\t", quote = F)

# Create correspondance table, i.e. which signature in which file corresponds to which consensus signatures
correspondanceTable <- data.frame(filename = gsub("(.+)\\|.+", "\\1", colnames(Wjoint)),
                                 originalSignature = as.numeric(gsub(".+\\|V(\\d+)", "\\1", colnames(Wjoint))) - 1,
                                 consensusSignature = heatmapResult$membership[colnames(Wjoint)],
                                 stringsAsFactors = F)
write.table(correspondanceTable, paste0(outPrefix, "_correspondanceSignatures.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

# And create a logical presence table of signature per tissue
presenceTable <- matrix(F, nrow = max(correspondanceTable$consensusSignature), ncol = length(unique(correspondanceTable$filename)))
colnames(presenceTable) <- unique(correspondanceTable$filename)
for(i in 1:nrow(correspondanceTable))
    presenceTable[as.numeric(correspondanceTable$consensusSignature[i]), correspondanceTable$filename[i]] <- T

colnames(presenceTable) <- gsub("_W", "", colnames(presenceTable))

write.table(presenceTable, paste0(outPrefix, "_presenceSignatures.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

# Plot histogram of presence absecence
presenceTable<- t(apply(presenceTable, 2, as.numeric))
pdf(paste0(outPrefix, "_presenceSignatures.pdf"), height = 9, width = 3)
heatmap.2(presenceTable, Rowv = T, Colv = T, col = c("grey", "#1b57b7"), trace = "none", key = F,
          colsep = 1:ncol(presenceTable), rowsep = 1:nrow(presenceTable))
dev.off()


#--------------------------------------------

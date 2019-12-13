# Creates a scatter plot of gene expression vs number of mutations
# Usage
# Rscript mutationVsExpression.R outHist.pdf mapFile1.txt [mapFile2.txt] [...]

source("../../R/mutationGeneAnnotation.R", chdir = T)
source("../../R/gtex.R", chdir = T)
source("../../R/plots.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
library("reshape")
library("dplyr")

cmdArgs <- commandArgs(T)
outplot <- cmdArgs[1]
mutationFiles <- cmdArgs[-1]

#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.7/"
#mutationFiles <- file.path(workingDir, list.files(workingDir))

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

appendCounts <- function(x, countFile, gtexId, colName) {

    geneExp <- readAllGtexExpression(gtexId, countFile)
    
    if(ncol(geneExp) < 2 | nrow(geneExp) < 2)
        return(NULL)
    
    geneExp$Name <- gsub("\\..+", "", geneExp$Name)
    rownames(geneExp) <- geneExp$Name
    geneExp <- geneExp[rownames(geneExp) %in% x$gene_id,]
    geneExp <- melt(geneExp)
    rownames(geneExp) <- paste0(geneExp[,1], ".", geneExp[,2])

    x$new <- geneExp[ paste0(x$gene_id, ".", x$gtexId), 3]
    x <- x[!is.na(x$new),]
    x <- x[x$new != 0,]
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

#' Correlates 

#' Creates histogram of correlations per sample including random sampling
histCors <- function(mutationsPerGene) {
    
    mutationsPerGene <- as_tibble(mutationsPerGene)
    mutationsPerGeneRandom <- randomizePerSample(mutationsPerGene, "mutations")
    
    mutationsPerGeneRandom$type <- "Random" 
    mutationsPerGene$type <- "Real" 
    
    mutationsPerGene <- rbind(mutationsPerGene, mutationsPerGeneRandom)
    mutationsPerGene$mutationsPerBp <- mutationsPerGene$mutations / geneLengths[mutationsPerGene$gene_id]
    
    toPlot <- 
        mutationsPerGene %>%
        group_by(sample, type) %>%
        summarise(cor = cor(mutationsPerBp, TPM))
        
    
    
    ggplot(toPlot, aes(x = cor, fill = type)) +
        geom_density(alpha = 0.3) + 
        geom_vline(aes(xintercept = median, colour = type), 
                   linetype = "dashed",
                   data = toPlot %>% group_by(type) %>% summarise(median = median(cor))) + 
        scale_fill_manual(values = c("#939dad", "#045be8")) +
        scale_colour_manual(values = c("#939dad", "#045be8")) +
        theme_bw()
}
    

#----------------------------------#
# MAIN
#----------------------------------#

####
# Read mutations and process mutations
mutationsAll <- list()
for(i in 1:length(mutationFiles)) {
    
    mutationFile <- mutationFiles[i]
    
    currentSample <- gsub(".txt", "", basename(mutationFile))
    mutations <- read.table(mutationFile, sep = "\t", stringsAsFactors = F)
    mutations$sample <- currentSample
    
    mutationsAll[[i]] <- mutations
    
}

mutationsAll <- do.call(rbind, mutationsAll)
mutationsAll <- annotateMutations(mutationsAll)

# Calculate number of mutations per gene
mutationsPerGene <- getMutationsPerGene(mutationsAll)

# Normalizing by gene length
geneLengths <- unique(mutationsAll[,c("gene_id", "cds_length")])
geneLengths <- setNames(geneLengths[,2], geneLengths[,1])
mutationsPerGene$mutationsPerBp <- mutationsPerGene$mutations / geneLengths[mutationsPerGene$gene_id]

# Adding gtexIds
gtexId <- setNames(sraToGtex(unique(mutationsPerGene$sample), formatOut = "long"), unique(mutationsPerGene$sample))
mutationsPerGene$gtexId <- gtexId[mutationsPerGene$sample]

#######
# Adding gene counts and TPM
mutationsPerGene <- appendCounts(mutationsPerGene,  CONFIG$auxiliaryFiles$expresionAllTissues, mutationsPerGene$gtexId, "TPM")
mutationsPerGene <- mutationsPerGene[ mutationsPerGene$mutations >= 3 & mutationsPerGene$TPM >= 1,]

if(is.null(mutationsPerGene)) {
    
    ggsave(outplot, ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100))
    
} else {


    mutationsPerGene <- appendCounts(mutationsPerGene,  CONFIG$auxiliaryFiles$readsGenesAllTissues, mutationsPerGene$gtexId, "read_counts")
    mutationsPerGene <- mutationsPerGene[mutationsPerGene$read_counts != 0,]

    ######
    # Plotting
    #p1 <- scatter(mutationsPerGene, alpha = 0.3, x = "TPM", y = "mutationsPerBp", regression = T, ylab = "Mutations per Bp", xlab = "TPM") + 
    #theme_noGrid()
    #p2 <- scatter(mutationsPerGene, alpha = 0.3, x = "read_counts", y = "mutations", regression = T, ylab = "Mutations", xlab = "Reads") + 
    #theme_noGrid()

    # Plotting histograms of correlation between expression and mutations,
    # it plots the random version too
    pHist <- histCors(mutationsPerGene)

    ggsave(outplot, pHist)
    #ggsave("~/plot.pdf", p1)
    #ggsave("~/plot2.pdf", p2)
}

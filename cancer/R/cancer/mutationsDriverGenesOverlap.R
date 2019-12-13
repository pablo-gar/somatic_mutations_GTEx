# Calculates the signifcance of the number of mutations in driver genes based
# on a comparison to a set of genes with similar expression patterns

# Usage
# Rscript mutationsDriverGenesOverlap.R out_ksTest.txt out_permutations.txt out_perm_byGene.txt out_perm_byGeneRandom.txt out_plotPerms.pdf /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/tcga_driverCancerGenes.txt /scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.5/*

source("../../../R/mutationGeneAnnotation.R", chdir = T)
source("../../../R/gtex.R", chdir = T)
source("../../../R/misc.R", chdir = T)
source("../../../R/geneTools.R", chdir = T)
library("reshape")
library("ggplot2")
library("rSubmitter")
miscSource <- "~/scripts/FraserLab/somaticMutationsProject/R/misc.R" # For rSubmitter

cmdArgs <- commandArgs(T)
output_ksTest <- cmdArgs[1]
output_permutations <- cmdArgs[2]
output_byGene <- cmdArgs[3]
output_byGeneRandom <- cmdArgs[4]
output_plotPermutations <- cmdArgs[5]
driverGeneFile <- cmdArgs[6]
mutationFiles <- cmdArgs[-(1:6)]


nPerm <- 1000

#workingDir <-  "/scratch/users/paedugar/somaticMutationsProject/mutationCount/map/Liver/n6_0.0_0.5/"
#mutationFiles <- file.path(workingDir, list.files(workingDir))
#driverGeneFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/tcga_driverCancerGenes.txt"

#----------------------------------#
# METHODS
#----------------------------------#
getMutationsPerGene <- function(x) {
    
    # Adding gtexIds
    gtexId <- setNames(sraToGtex(unique(x$sample), formatOut = "long"), unique(x$sample))
    x$gtexId <- gtexId[x$sample]
    
    result <- by(x, x$gtexId, function(x) {
                 y <-as.data.frame(table(x$gene_id), stringsAsFactors = F)
                 colnames(y) <- c("gene_id", "mutations")
                 y$gtexId <- x$gtexId[1]
                 y$sample <- x$sample[1]
                 rownames(y) <- y$gene_id
                 return(y)
                })
    return(do.call(rbind,result))
}

#' Reads the expression of a vector of genes
readExpresison <- function(countFile, gtexId) {

    geneExp <- readAllGtexExpression(gtexId, countFile)
    
    if(ncol(geneExp) < 2 | nrow(geneExp) < 2)
        return(NULL)
    
    geneExp$Name <- gsub("\\..+", "", geneExp$Name)
    rownames(geneExp) <- geneExp$Name

    return(geneExp)
    
}

#' Merges expression and mutations data.frames
#' it requires tha both table have equivalent rownames
mergeExpressionMutation <- function(genes, geneExp, mutations) {
    
    # Melts gene expression only with genes of interest
    geneExp <- geneExp[rownames(geneExp) %in% genes,]
    geneExp <- melt(geneExp)
    rownames(geneExp) <- paste0(geneExp[,2], ".", geneExp[,1])
    colnames(geneExp) <- c("gene_id", "gtexId", "expression")
    
    # Merges with mutation
    geneExp$mutations <- mutations[rownames(geneExp), "mutations"]
    geneExp[is.na(geneExp$mutations), "mutations"] <- 0

    return(geneExp)
}

#' Gets number of mutations of a set of genes
#' it requires tha both table have equivalent rownames
#' only those present in the geneExp table
getMutationNumber <- function(genes, geneExp, mutations) {
    
    # Melts gene expression only with genes of interest
    genes <- genes[genes %in% rownames(geneExp)]
    
    # Creating vector of all combinations of individual/genes
    geneInd <- expand.grid(unique(mutations$gtexId), unique(genes))
    
    # Getting mutations for the genes
    geneInd$mutations <-  mutations[paste0(geneInd[,1], ".", geneInd[,2]), "mutations"]
    colnames(geneInd) <- c("gtexId", "gene_id", "mutations")
    geneInd$mutations[is.na(geneInd$mutations)] <- 0
    
   # result <- mutations[paste0(geneInd[,1], ".", geneInd[,2]), "mutations"]
   # result[is.na(result)] <- 0
   # return(result)
    
    return(geneInd)
}


#' Given a vector of genes and a expression matrix of n_genes x n_individuals
#' returns a vector of genes that fall within the same corresponding quantiles
#' for each original gene
selectSimilarExpression <- function(genes, geneExp, nQuant = 100) {
    
    avgExp <- rowMeans(geneExp[, !colnames(geneExp) %in% "Name"])
    avgExp <- avgExp[ avgExp != 0 ]
    quants <- quantile(avgExp, probs = seq(0, 1, length.out = nQuant + 1))
    allGenesIntervals <- findInterval(avgExp, quants)
    
    intervals <- findInterval(avgExp[genes], quants)
    intervals[is.na(intervals)] <- min(allGenesIntervals)
    for(i in 1:length(intervals)) {
        randomI <- which(intervals[i] == allGenesIntervals)
        newGene <- genes[i] 
        while(newGene %in% genes) {
            newGene <- names(avgExp[sample(randomI, 1)])
        }
        genes[i] <- newGene
    }
    
    return(genes)
    
}

singlePermutation <- function(x, mutationsDriver, mutationsDriverRandom, uniqDrivers, randomDrivers,  geneExp, mutationsPerGene,  newRandomPerIteration = F) {
                                  
    printTime(paste0("Permutation ", x, "\n"))
    
    pvals <- data.frame(pvalReal = 1, pvalRandom = 1, realMeanDiff = 0, randomMeanDiff = 0)

    randomGenes <- selectSimilarExpression(uniqDrivers, geneExp)
    
    if(newRandomPerIteration)
        mutationsDriverRandom <- getMutationNumber(selectSimilarExpression(uniqDrivers, geneExp), geneExp, mutationsPerGene)
    
    # Getting number of mutations for random groups
    mutationsRandom <- getMutationNumber(randomGenes, geneExp, mutationsPerGene)
    
    # Performing Tests
    pvals$pvalReal[1] <- wilcox.test(mutationsDriver$mutations, mutationsRandom$mutations)[["p.value"]]
    pvals$pvalRandom[1] <- wilcox.test(mutationsDriverRandom$mutations, mutationsRandom$mutations)[["p.value"]]
    
    pvals$realMeanDiff[1] <- log2(mean(mutationsDriver$mutations) / mean(mutationsRandom$mutations))
    pvals$randomMeanDiff[1] <- log2(mean(mutationsDriverRandom$mutations) / mean(mutationsRandom$mutations))
    
    # Doing gene-by-gene analysis
    genesDiff <- log2(
                      (tapply(mutationsDriver$mutations, mutationsDriver$gene_id, function(x) sum(x)) + 0.001) / 
                      (tapply(mutationsRandom$mutations, mutationsRandom$gene_id, function(x) sum(x)) + 0.001)
                     )
    
    genesDiffRandom <- log2(
                            (tapply(mutationsDriverRandom$mutations, mutationsDriverRandom$gene_id, function(x) sum(x)) + 0.001) / 
                            (tapply(mutationsRandom$mutations, mutationsRandom$gene_id, function(x) sum(x)) + 0.001)
                           )
    
    
    return(list(pvals, genesDiff, genesDiffRandom))
}
    


#----------------------------------#
# MAIN
#----------------------------------#

###
# Read cancer genes
driverGenes <- read.delim(driverGeneFile, sep = "\t", stringsAsFactors = F, header = T)
driverGenes$gene_id <- aliasToEnsembl(driverGenes$Gene)

####
# Read mutations and process mutations
mutationsAll <- list()
for(i in 1:length(mutationFiles)) {
#for(i in 1:3) {
    
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

#######
# Read expression 
geneExp <- readExpresison(CONFIG$auxiliaryFiles$readsGenesAllTissues, mutationsPerGene$gtexId)

####
# Select only data for which gene drive info is available

# select drivers with expression info
driverGenes <- driverGenes[ driverGenes$gene_id %in% geneExp$Name,]
uniqDrivers <- unique(driverGenes$gene_id)


####
# PERFORM TESTS
# 1. Performs wilcox test of driver gene mutations vs randomly selected genes (same expression levels)
#   - Repeats 1000 times selecting new set of genes each time
# 2. Performs wilcox test of randomly selected gene mutations vs randomly selected genes (same expression levels)
#   - Repeats 1000 times selecting new set of genes each time (only second list)

randomDrivers <- selectSimilarExpression(uniqDrivers, geneExp)

# Getting number of mutations
mutationsDriver <-  mergeExpressionMutation(uniqDrivers, geneExp, mutationsPerGene)
mutationsDriverRandom <- mergeExpressionMutation(randomDrivers, geneExp, mutationsPerGene)

# Doing permutations
permResults <- "failed"
counter <- 1
while(is.character(permResults)) {
    if (counter > 4)
        stop("Tried 4 times running supper apply")
    flush.console()
    cat("Trying ", counter, " try")
    permResults <- tryCatch(superApply(1:nPerm, singlePermutation,
                                       # Super apply pars
                                       tasks = 10, partition = "hbfraser,owners,hns", time = "3:00:00", mem = "8G", proc = 1, 
                                       packages = "", clean = T,
                                       workingDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles", 
                                       # Function pars
                                       mutationsDriver = mutationsDriver, mutationsDriverRandom = mutationsDriverRandom, 
                                       uniqDrivers = uniqDrivers, randomDrivers = randomDrivers, newRandomPerIteration = T, 
                                       geneExp = geneExp, mutationsPerGene = mutationsPerGene),
                            error = function(err) {
                                "failed"
                            })
    counter <- counter + 1
}



pvals <- do.call(rbind, lapply(permResults, function(x) x[[1]]))
byGene <- do.call(rbind, lapply(permResults, function(x) x[[2]]))
byGeneRandom <- do.call(rbind, lapply(permResults, function(x) x[[3]]))
    
# Comparing distributions of pvalues with kolmogorov test
results <- data.frame(nPerm = nPerm, 
                      ksPvalue = ks.test(pvals$real, pvals$random)[["p.value"]],
                      wilcoxPairedPvalue = wilcox.test(pvals$real, pvals$random, paired = T, alternative = "less")[["p.value"]]
                      )

# Comparing distributions of pvalues with paired wilcox test

#Plotting distributions
p <- ggplot(melt(pvals[,c("pvalRandom", "pvalReal")]), aes(x = value, fill = variable)) +
         geom_density(alpha = 0.5) +
         theme_bw()
     
# Saving outputs
write.table(results, output_ksTest, sep = "\t", quote = F, row.names = F)
write.table(pvals, output_permutations, sep = "\t", quote = F, row.names = F)
write.table(byGene, output_byGene, sep = "\t", quote = F, row.names = F)
write.table(byGeneRandom, output_byGeneRandom, sep = "\t", quote = F, row.names = F)
ggsave(output_plotPermutations, p) 

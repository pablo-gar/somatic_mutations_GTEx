##------------------------------------------------------------
## INFO
## Analyze the pairwise pileup comparisons
##
## Author: Pablo E. Garcia-Nieto
##
## Date: 08/16/2017
##
## Example
##
## Rscript pileupComparisonsAnalysis.R ~/data/bloodGTExProject/pileup_pairwise_comparison/test_MAF_1 /scratch/PI/hbfraser/gtex/anno/GtexSraRunTable.txt ~/results/FraserLab/skinGTExProject/bloodMethodValidation/generalCounts/today/test_MAF_1
##------------------------------------------------------------


##------------------------------------------------------------
## PARAMETERS

args <- commandArgs(TRUE)
pileupDir <- args[1]  #"~/data/bloodGTExProject/pileup_pairwise_comparison"
outfile <- args[2]  #"~/results/FraserLab/skinGTExProject/bloodMethodValidation/generalCounts/today"

#pileupDir <-"~/data/bloodGTExProject/SNV/MAF_1/RNASEQ"
#outdir <- "~/results/FraserLab/skinGTExProject/bloodMethodValidation/generalCounts/today/MAF_1/SNVs/RNASEQ"

# Does nSupporting cutoff have to be met in both files?
nSupporting_both <- F

##------------------------------------------------------------


##------------------------------------------------------------
## GLOBAL

MUT_KEY <- matrix(ncol = 2, byrow = T,
                  data = c(
                    "C.T", "G.A",
                    "C.A", "G.T",
                    "C.G", "G.C",
                    "T.A", "A.T",
                    "T.C", "A.G",
                    "T.G", "A.C"
                  )
                )

# Indeces for reference and alternate alleles counts, for each experiment 
REF_I <- 6
ALT_I <- 7

# Indeces for reference and alternate allele names
REF_A <- 3
ALT_A <- 4

# PLOT OPTIONS
COL_1 <- "#C8C8A9"
COL_2 <- "#F1C7A8"

WIDE_PLOT <- list(width = 10, height = 5)
LARGE_PLOT <- list(width = 10, height = 10)

##------------------------------------------------------------



##------------------------------------------------------------
## METHODS

getTableFile <- function(x) {
    
    # Returns a data frame where each columns has info for 
    # each pairwsie-pileup comparison file found in the pileupDir
    #
    # @args x - string - path to location of pileup files
    #
    # @return - data.frame - containing the info of each pileup file:
    #     - individual id
    #     - exp id
    #     - coverage cutoff
    #     - nSupporting reads cutoff
    #     - path to file
    
    pileupFiles <- list.files(pileupDir)
    pileupFiles <- pileupFiles[grep(".txt$", pileupFiles)]
    
    results <- data.frame(ind = "", sampleA = "", coverage = 0, nSupporting = 0, location = file.path(x, pileupFiles), stringsAsFactors = F)
                          
    pileupFiles <- gsub(".txt", "", pileupFiles)
    pileupFiles <- strsplit(pileupFiles, "_")
    
    for (i in 1:length(pileupFiles)) {
        results[i, "ind"] <- pileupFiles[[i]][1]
        results[i, "sampleA"] <- pileupFiles[[i]][1]#pileupFiles[[i]][1]
        results[i, "coverage"] <- 40 #pileupFiles[[i]][3]
        results[i, "nSupporting"] <- 4#pileupFiles[[i]][5]
    }
    
    return(results)
}

getPileups <- function(tableFile) {
    
    # Reads all pileups from the tableFile into a list of data.frames
    # appends to the pileups the coverage and nSupporting cutoffs as well as
    # individual id and experiments ids

    allPileups <- list()
    i <- 1
    for(coverage in unique(tableFile$coverage)) {
        for(nSupporting in unique(tableFile$nSupporting)){
            
            currentFiles <- tableFile[ tableFile$coverage == coverage & tableFile$nSupporting == nSupporting,]
            
            for(j in 1:nrow(currentFiles)) {
                currentFile <- currentFiles[j, "location"]
                
                if(file.size(currentFile) == 0) 
                    next
                
                currentPileup <- readPileup(currentFile)
                
                # Appending info
                currentPileup$coverage <- coverage
                currentPileup$nSupporting <- nSupporting
                currentPileup$ind <- currentFiles[j, "ind"]
                currentPileup$sampleA <- currentFiles[j, "sampleA"]
                
                allPileups[[i]] <- currentPileup
                
                i <- i + 1
            }
        }
    }
    
    return(allPileups)
}


readPileup <- function(x) {
    
    x <- read.table(x, sep = "\t", header = F, stringsAsFactors = F)
    
    return(x)
    
}

simplifyMutationTypes <- function(mutCount) {
    
    # From a data.frame of counts of SNVs, reduces the number or rows
    # by merging equivalent mutation type counts into one (e.g. C>T and G>A)
    # 
    # @args mutCount - data.frame - table of mutations counts, output of getCountsByType()
    # 
    # @return - data.frame - table of mutations counts with less row number
    
    results <- list()
    for (i in 1:nrow(MUT_KEY)) { 
        mutPair <- MUT_KEY[i,]
        currentMutCount <- mutCount[mutCount$mutType  %in% mutPair, , drop = F]
        currentMutCount <- currentMutCount[, -ncol(currentMutCount)]
        currentMutCount <- colSums(currentMutCount)
        results[[i]] <- currentMutCount
    }
    
    results <- as.data.frame(do.call(rbind, results))
    results$mutType <- MUT_KEY[,1]
    
    return(results)
}
        
    

appendTissue <- function(x, idColumns, sraTableFile) {
    
    sraTable <- read.table(sraTableFile, sep = "\t", stringsAsFactors = F, header = T)
    sraTable$Sample_Name_s <- gsub("_rep\\d+", "", sraTable$Sample_Name_s)
    
    sraTable <- sraTable[!duplicated(sraTable$Sample_Name_s),]
    
    
    rownames(sraTable) <- sraTable$Sample_Name_s
    
    for (col in idColumns) {
        x$extra <- sraTable[ x[,col], "body_site_s"]
        colnames(x)[ncol(x)] <- paste0(col, "Tissue")
    }
    
    return(x)
    
}
       
createHistogram <- function(pileupList) {
    
    require(ggplot2)
    
    ref_i <- REF_I
    alt_i <- ALT_I
    
    currentExp <- mergeDataExp(pileupList, ref_i, alt_i)
    #currentExp$MAF <- currentExp[ ,alt_i] /  (currentExp[ ,alt_i] + currentExp[ ,ref_i])
    currentExp$MAF <- currentExp[ ,alt_i] /  currentExp[ ,ref_i]
    
    
    currentExp$coverage <- factor(currentExp$coverage, levels = sort(as.numeric(as.character(unique(currentExp$coverage))), decreasing = F), ordered = T)
    currentExp$nSupporting <- factor(currentExp$nSupporting, levels = sort(as.numeric(as.character(unique(currentExp$nSupporting))), decreasing = F), ordered = T)
    
    legend <- data.frame( x = Inf, y = Inf,
                         label = paste0("Mean = ", mean(currentExp$MAF), "\n",
                                        "Min = ", signif(min(currentExp$MAF), 4), "\n",
                                        "Max = ", signif(max(currentExp$MAF), 4), "\n", 
                                        "Nind = ", length(unique(currentExp$ind)), "\n",
                                        "Nmut = ", nrow(currentExp)
                                        )
                        )
    
    p <- ggplot(currentExp, aes(x = MAF)) + 
    geom_density() + 
    coord_cartesian(xlim = c(0, 1)) + 
    facet_grid(nSupporting ~ coverage) +
    geom_text(aes(label = label, x = x, y = x), data = legend, vjust = 1, hjust = 1) + 
    theme_bw()  
        
    
    return(p)
    
}
        
mergeDataExp <- function(pileupList, ref_i, alt_i) {
    
    merged <- list()
    for (i in 1:length(pileupList)) {
        currentPileup <- pileupList[[i]]
        merged[[i]] <- currentPileup[ currentPileup[,alt_i] >= currentPileup[1, "nSupporting"], ]
    }
    
    merged <- do.call(rbind, merged)
    
    return(merged)
}

savePlots <- function(plotList, outfile, ggplotArgs = list()) {
    
    for (i in 1:length(plotList)) {
        
        #plotName <- names(plotList)[i]
        filename <- outfile #file.path(outdir, paste0(prefix, plotName, ".pdf"))
        
        args <- c(list(filename = filename, plot = plotList[[i]]), ggplotArgs)
        
        do.call(ggsave, args)
    }
    return()
}
        
    
    
##------------------------------------------------------------



##------------------------------------------------------------
## MAIN

tableFile <- getTableFile(pileupDir)
allPileups <- getPileups(tableFile)

# Creates MAF histogram
histAll <- list()
histAll[["SNV"]] <- createHistogram(allPileups)

# Saves histograms
savePlots(histAll, outfile, WIDE_PLOT)

##------------------------------------------------------------

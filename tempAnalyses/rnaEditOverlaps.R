# THIS HAS TO BE CALLED FROM WHERE IT RESIDES
#
# Takes all mutation maps from tissues, finds the overlaps with RNA edits and performs:
#
# -Percentages of RNA edits in samples
# -Type of mutations in the overalps
# -Type of mutations without the overalps

# Usage Rscript

library("rSubmitter")
library("purrr")
library("ggplot2")
source("../R/ggthemes.R")

# GLOBALS
BAR_WIDTH = 0.3

# PARS
rnaEditFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/Human_AG_all_hg19_v2.txt"

rootFolder <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount"
mapFolder <- file.path(rootFolder, "map")
countFolder <- file.path(rootFolder, "count")

maf <- "n6_0.0_0.5"

outDir <- "/scratch/users/paedugar/somaticMutationsProject/rnaEdit"
mapOutDir <- file.path(outDir, "mutationCount/map")

MUT_TO_BED <- "$HOME/scripts/FraserLab/somaticMutationsProject/bin/mutationToBed"
readCountFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/readCounts.txt_2018-05-14_18:26:06"

# Perform the overlap with RNA edits?
doRNAeditOverlap <- F

# Outdirs for analyses
generalAnalysesPrefix <- "generalMutationAnalyses/results"
percentageVariantsInRNAsites <- mutationTypesNormalized_outDir <- file.path(outDir, generalAnalysesPrefix, "percentageInRNA_editSites")
mutationTypesNormalized_outDir <- file.path(outDir, generalAnalysesPrefix, "countsPerMutTypeNormalized")
mutationTypes_outDir <- file.path(outDir, generalAnalysesPrefix, "countsPerMutType")
totalCountsBox_outDir <- file.path(outDir, generalAnalysesPrefix, "totalCountsBox")

#-------------------------------
# METHODS

#' Intersect a mutation map table and the RNA edit table
#' @param rnaEditFileBED - file in BED formatm needs to have one extra column from original (end coordinate)
intersectMut <- function(mutationFile, rnaEditFileBED, outfile) {

    MUT_TO_BED <- "$HOME/scripts/FraserLab/somaticMutationsProject/bin/mutationToBed"
    
    tempMutBed <- paste0(outfile, ".temp")

    # Converts mutation File to bed
    system(paste(MUT_TO_BED, mutationFile, ">", tempMutBed), intern = T)

    # Constructs bedtools intersect
    cmd <- paste("bedtools intersect -a", tempMutBed, "-b", rnaEditFileBED, ">", outfile)
    system(cmd, intern = T)

    file.remove(tempMutBed)

    return(invisible(NULL))
}

#' Eliminates third column of a file
eliminate3rdColumn <- function(x) {

    firstFile <- paste0(x, ".temp1")
    secondFile <- paste0(x, ".temp2")

    system(paste("cut -f 1-2", x, ">", firstFile))
    system(paste("cut -f 4-", x, ">", secondFile))
    system(paste("paste", firstFile, secondFile, "> temp && mv temp", x))


    file.remove(firstFile)
    file.remove(secondFile)
    return(invisible(NULL)) 
}
     
#' Print mutations in RNA edit sites for one mut table
getMutationsInRnaEdit <- function(x, tissueOutDir, tempRnaEditBed) {

    mutTableOut <- file.path(tissueOutDir, basename(x))

    intersectMut(x, tempRnaEditBed, mutTableOut)
    eliminate3rdColumn(mutTableOut)

    return(invisible(NULL))
}

#' Print mutations in RNA edit sites for all mut table in a dir using parallel computing
getMutationsInRnaEditFolder <- function(tissueDir, tempRnaEditBed, mapOutDir) {
    
    tissueOutDir <- file.path(mapOutDir, basename(dirname(tissueDir)))
    dir.create(tissueOutDir,showWarnings = F, recursive = T)
    
    
    # Working only with samples that do not exist or whose files are empty
    outAllMutTables <- file.path(tissueOutDir, list.files(tissueOutDir))
    info <- file.info(outAllMutTables)
    nonEmpty <- rownames(info[info$size > 0, ])
    nonEmpty <- basename(nonEmpty)
    
    # Gathering the input files
    inputSamples <- list.files(tissueDir)
    inputSamples <- inputSamples[!inputSamples %in% nonEmpty]
    allMutTables <-  file.path(tissueDir, inputSamples)
    
    
    if (length(allMutTables) != 0) {
        superApply(allMutTables, getMutationsInRnaEdit, tissueOutDir = tissueOutDir, tempRnaEditBed = tempRnaEditBed, 
                   partition = "hbfraser,owners,hns,normal", mem = "4G", time = "2:00:00", proc = 1,
                   tasks = length(allMutTables), workingDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles",
                   packages = "")
    }
    
    return(invisible(NULL))
    
}

#' Print and count mutation types
#' takes fileVector, where the first element is the input map file and the second is the output count file
printCounts <- function(fileVector) {
    
    mapFile <- fileVector[1]
    outCountFile <- fileVector[2]
    
    dir.create(dirname(outCountFile), showWarnings = F, recursive = T)
    
    outMat <- matrix(0, nrow = 2, ncol = 4)
    
    flush.console()
    cat(file.info(mapFile)$size, " ", mapFile, "\n")
    
    mapTable <- tryCatch(read.delim(mapFile, sep = "\t", stringsAsFactors = F), error = function(cond) data.frame())
    
    if(nrow(mapTable) != 0) {
    
    
        frequencies <- table(paste0(mapTable[,3], mapTable[,4]))
        
        outMat[1,1] <- sum(frequencies[names(frequencies) %in% c("TA", "AT")])
        outMat[1,3] <- sum(frequencies[names(frequencies) %in% c("TG", "AC")])
        outMat[1,4] <- sum(frequencies[names(frequencies) %in% c("TC", "AG")])
        outMat[2,1] <- sum(frequencies[names(frequencies) %in% c("CA", "GT")])
        outMat[2,2] <- sum(frequencies[names(frequencies) %in% c("CT", "GA")])
        outMat[2,3] <- sum(frequencies[names(frequencies) %in% c("CG", "GC")])
        
        write.table(outMat, outCountFile, sep = "\t", row.names = F, col.names = F, quote = F)
        return(outCountFile)
    }
    
    
    return(NULL)
}
    
    
    



#-------------------------------


#-------------------------------
# MAIN




##########
# FINDING RNA EDIT OVERLAPS AND PRINTING MAPS

allTissues <- file.path(list.dirs(mapFolder, recursive = F), maf)
    
if(doRNAeditOverlap) {
    dir.create(mapOutDir, showWarnings = F, recursive = T)
    
    # Converts rna edit file to bed
    tempRnaEditBed <- paste0(rnaEditFile, ".temp")
    system(paste(MUT_TO_BED, rnaEditFile, "| tail -n +2 >", tempRnaEditBed))

        
    nullList <- superApply(allTissues, getMutationsInRnaEditFolder, tempRnaEditBed = tempRnaEditBed, mapOutDir = mapOutDir,
                           partition = "hbfraser,owners,hns,normal", mem = "4G", time = "4:00:00", proc = 1,
                           tasks = length(allTissues), workingDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles",
                           packages = "rSubmitter")
    
    file.remove(tempRnaEditBed)
}


#######
# CREATING COUNT MATRICES
allMapFiles <- unlist(map(list.dirs(mapOutDir, recursive = F), ~ file.path(.x, list.files(.x)) ))
allCountFiles <- gsub("/map/", "/count/", allMapFiles)

fileList <- map2(allMapFiles, allCountFiles, ~ c(.x, .y))

allCountFiles <- superApply(fileList, printCounts, 
                       partition = "hbfraser,owners,hns,normal", mem = "2G", time = "10:00", proc = 1,
                       tasks = 100,  workingDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles",
                       packages = "")


allCountFiles <- unlist(allCountFiles)
allCountFiles <- allCountFiles[!is.null(allCountFiles)]



#######
# DOING ANALYSES
# Here just to play around with rSubmmitter I'm going to do individual submissions

##
# Mutation types
allCountTissues <- unique(dirname(allCountFiles))
dir.create(mutationTypes_outDir, showWarnings = F, recursive = T)

# Submits jobs
allJobs <- map(allCountTissues, function(x, mutationTypes_outDir) {
                        tissue <- basename(x)
                        outplot <- file.path(mutationTypes_outDir, paste0(tissue, ".pdf"))
                        cmd <- paste0("module load R/3.4.0; cd ~/scripts/FraserLab/somaticMutationsProject/generalMutationAnalyses/R/; Rscript mutationTypesBars.R ", outplot, " ", x, "/*")
                        job <- Job$new(cmd, partition= "hbfraser,owners,hns,normal", mem = "2G", time = "10:00", proc = 1,
                                       outDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles")
                        job$submit()
                        return(job)
                       }, mutationTypes_outDir = mutationTypes_outDir)

# Waits for jobs
walk(allJobs, ~ .x$wait())

# Cleans for jobs
walk(allJobs, ~ .x$clean())


##
# Mutation types normalized
dir.create(mutationTypesNormalized_outDir, showWarnings = F, recursive = T)

allJobs <- map(allCountTissues, function(x, mutationTypes_outDir) {
                        tissue <- basename(x)
                        outplot <- file.path(mutationTypes_outDir, paste0(tissue, ".pdf"))
                        cmd <- paste0("module load R/3.4.0; cd ~/scripts/FraserLab/somaticMutationsProject/generalMutationAnalyses/R/; Rscript mutationTypesBarsNormilized.R ", outplot, " ", x, "/*")
                        job <- Job$new(cmd, partition= "hbfraser,owners,hns,normal", mem = "2G", time = "10:00", proc = 1,
                                       outDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles")
                        job$submit()
                        return(job)
                       }, mutationTypes_outDir = mutationTypesNormalized_outDir)

# Waits for jobs
walk(allJobs, ~ .x$wait())
# Cleans for jobs
walk(allJobs, ~ .x$clean())


##
# Boxplots normilized
dir.create(totalCountsBox_outDir, showWarnings = F, recursive = T)
allCountTissuesRoot <- unique(dirname(allCountTissues))

# This is a one-step loop, daah. But I kept it to encapsulate variables and becasue I'm lazy to change it
allJobs <- map(allCountTissuesRoot, function(x, totalCountsBox_outDir) {
               
                        # Normalized by read counts
                        outplot <- file.path(totalCountsBox_outDir, "normalized.pdf")
                        cmd <- paste0("module load R/3.4.0; cd ~/scripts/FraserLab/somaticMutationsProject/generalMutationAnalyses/R/; Rscript mutationCountBoxplotTissuesNormalized.R ", 
                                      x, " '' ", readCountFile, " ", outplot)
                        
                        # Raw
                        outplot <- file.path(totalCountsBox_outDir, "raw.pdf")
                        cmd2 <- paste0("Rscript mutationCountBoxplotTissues.R ", 
                                       x, " '' ", outplot)
                        
                        cmd <- c(cmd, cmd2)
                        
                        job <- Job$new(cmd, partition= "hbfraser", mem = "8G", time = "20:00", proc = 1,
                                       outDir = "/scratch/users/paedugar/somaticMutationsProject/clusterFiles")
                        job$submit()
                        return(job)
                       }, totalCountsBox_outDir = totalCountsBox_outDir)

# Waits for the single job
walk(allJobs, ~ .x$wait())

# Cleans for jobs
walk(allJobs, ~ .x$clean())


##
# Calculate percentages of variants called in RNA edit sites

percentages <- map_dbl(allCountTissues, function(tissueFolder, originalCountFolder, maf) {
                       
                  countFiles <- file.path(tissueFolder, list.files(tissueFolder))
                  names(countFiles) <- basename(countFiles)
                  countsAll <- map_dfr(countFiles, function(countFile, originalCountFolder, maf) {
                                        
                                        sampleId <- basename(countFile)
                                        tissue <- basename(dirname(countFile))
                                        originalCountFile <- file.path(originalCountFolder, tissue, maf, sampleId)
                                        
                                        rnaEdits <- as.matrix(read.table(countFile, sep = "\t", header = F))
                                        allVariants <- as.matrix(read.table(originalCountFile, sep = "\t", header = F))
                                        
                                        return (data.frame(rnaEdits = sum(rnaEdits), allVars = sum(allVariants)))
                                        }, originalCountFolder = originalCountFolder, maf = maf
                                        )
                  
                  return(sum(countsAll[,1]) / sum(countsAll[,2]) * 100)
                  }, originalCountFolder = countFolder, maf = maf
                  )

toPlot <- data.frame(Percent_variants_called_in_RNA_edit_sites = percentages, tissue = basename(allCountTissues))
toPlot$tissue <- factor(toPlot$tissue, ordered = T, levels = toPlot$tissue[order(toPlot[,1])])

p <- ggplot(toPlot, aes(x = tissue, y = Percent_variants_called_in_RNA_edit_sites)) + 
geom_bar(stat = "identity") +
coord_flip() +
ylim (0, max(toPlot$Percent_variants_called_in_RNA_edit_sites + 20)) + 
theme_grid_x()


dir.create(percentageVariantsInRNAsites, showWarnings = F, recursive = T)
ggsave(file.path(percentageVariantsInRNAsites, "allTissues.pdf"), p, height = nrow(toPlot) * BAR_WIDTH)




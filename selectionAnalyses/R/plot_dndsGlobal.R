library(ggplot2)
source("../../R/ggthemes.R", chdir = T)
source("../../R/plots.R", chdir = T)

args <- commandArgs(T)

workingDir <- args[1] # Dir with selection results, e.g "/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/"
maf <- args[2] # Maf of interest 
plotPrefix <- args[3]

#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/"
#maf <- "n6_0.0_0.5"
#plotPrefix <- "/scratch/users/paedugar/somaticMutationsProject/selection/dndsGlobalPlots/dnds_"

WIDTH_PER_FACET <- 6 # In inches
HEIGHT <- 4.5

if(!dir.exists(workingDir))
    stop(paste(workingDir, " does not exist"))


tissues <- basename(list.dirs(workingDir, recursive = F)) 
if(length(tissues) == 0)
    stop(paste(workingDir, " has no tissue folders"))

# Reading data
dndsoutAll <- list()
for(tissue in tissues) {
        
    tissuePath <- file.path(workingDir, tissue, maf)
    if(!dir.exists(tissuePath))
        next
    
    if(grepl("EXO", tissuePath))
        next
        
    globFiles = list.files(tissuePath)
    globFiles = globFiles[grep("glob", globFiles)]
    if(length(tissues) == 0)
        stop(paste(tissuePat, " has no global selection files"))
        
        
    for (i in 1:length(globFiles)) {
        samplePath <- file.path(tissuePath, globFiles[i])
        dndsout <- read.table(samplePath, sep = "\t", stringsAsFactors = F, header = T)
        dndsout$tissue <- tissue
        dndsout$expression <- gsub(".+,geneExpression:(.+?),.+", "\\1", globFiles[i])
        dndsout$maf <- gsub(".+,selectionMaf:(.+?),.+", "\\1", globFiles[i])
        dndsoutAll[[paste0(i, tissue)]] <- dndsout
    }   
        
}

dndsoutAll <- do.call(rbind, dndsoutAll)
dndsoutAll$cihigh[dndsoutAll$cihigh > 4] <- 0


# Comparing maf
allPlots <- list()
pointRangePars <- list(x = "tissue", y = "mle", errorBarYmin = "cilow", errorBarYmax = "cihigh", pSize = 0.4)
facet_x <- "maf"
groupBy <- "expression"

for(selType in unique(dndsoutAll$name)) {
    toPlotGlobal <-  dndsoutAll[dndsoutAll$name == selType,]
    allPlots[[selType]] <- list()
    
    for(secondType in unique(dndsoutAll[ ,groupBy])) {
        
        toPlot <- toPlotGlobal[toPlotGlobal[ ,groupBy] == secondType,]
        
        # Order tissues based on average
        tissuesToOrder <- toPlot[ toPlot[,facet_x] == "0%_100%", c("mle", "tissue")]
        toPlot$tissue <- factor(toPlot$tissue, levels = tissuesToOrder[order(tissuesToOrder$mle), "tissue"])
    
        p <- do.call(pointRange, c(pointRangePars, list(dataframe = toPlot, facet_x = facet_x, scales = "free_x")))
        p <- p + theme_grid_y() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
        
        # Adding median line per facet
        medianLine <- tapply(toPlot$mle, toPlot[,facet_x], median)
        medianLine <- data.frame(yintercept = medianLine, facet_x = names(medianLine))
        colnames(medianLine)[ncol(medianLine)] <- facet_x
        
        p <- p + geom_hline(aes(yintercept = yintercept), data = medianLine, linetype = "dashed", colour = "lightskyblue3")
        
        allPlots[[selType]][[secondType]] <- p
        
        
        filename <- paste0(plotPrefix, "dotPlot_", selType, "_", groupBy, "_", secondType, "_facet_", facet_x, ".pdf")
        filename <- gsub("\\%", "", filename)
        flush.console()
        ggsave(filename, allPlots[[selType]][[secondType]], height = HEIGHT, width = WIDTH_PER_FACET * length(unique(toPlot[,facet_x])))
    }
}


allPlots <- list()
pointRangePars <- list(x = "tissue", y = "mle", errorBarYmin = "cilow", errorBarYmax = "cihigh", pSize = 0.4)
facet_x <- "expression" 
groupBy <- "maf"

for(selType in unique(dndsoutAll$name)) {
    toPlotGlobal <-  dndsoutAll[dndsoutAll$name == selType,]
    allPlots[[selType]] <- list()
    
    for(secondType in unique(dndsoutAll[ ,groupBy])) {
        
        
        toPlot <- toPlotGlobal[toPlotGlobal[ ,groupBy] == secondType,]
        
        # Order tissues based on average
        tissuesToOrder <- toPlot[ toPlot[,facet_x] == "AllExpressed", c("mle", "tissue")]
        toPlot$tissue <- factor(toPlot$tissue, levels = tissuesToOrder[order(tissuesToOrder$mle), "tissue"])
    
        p <- do.call(pointRange, c(pointRangePars, list(dataframe = toPlot, facet_x = facet_x, scales = "free_x")))
        p <- p + theme_grid_y() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
        
        # Adding median line per facet
        medianLine <- tapply(toPlot$mle, toPlot[,facet_x], median)
        medianLine <- data.frame(yintercept = medianLine, facet_x = names(medianLine))
        colnames(medianLine)[ncol(medianLine)] <- facet_x
        
        p <- p + geom_hline(aes(yintercept = yintercept), data = medianLine, linetype = "dashed", colour = "lightskyblue3")
        
        allPlots[[selType]][[secondType]] <- p
        
        filename <- paste0(plotPrefix, "dotPlot_", selType, "_", groupBy, "_", secondType, "_facet_", facet_x, ".pdf")
        filename <- gsub("\\%", "", filename)
        ggsave(filename, allPlots[[selType]][[secondType]], height = HEIGHT, width = WIDTH_PER_FACET * length(unique(toPlot[,facet_x])))
        
    }
}



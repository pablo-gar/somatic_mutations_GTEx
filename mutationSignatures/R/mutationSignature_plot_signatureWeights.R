if(!dir.exists("../../R"))
    stop("Can't find the shared R path, make sure to run script from directory where it lives")

library(gplots)
source("../../R/gtex.R", chdir = T)

# Gathering arguments
args <- commandArgs(T)
Hfile <- args[1] #"/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Liver-n6_0.0_0.5_H.txt"
outFile <- args[2] #"/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Liver-n6_0.0_0.5_H_plot.pdf"

#Hfile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Liver-n6_0.0_0.5_H.txt"
#outFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Liver-n6_0.0_0.5_H_plot.pdf"
#Hfile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Skin_Sun_Exposed_Lower_leg-n6_0.0_0.5_H.txt"
#outFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Skin_Sun_Exposed_Lower_leg-n6_0.0_0.5_H_plot.pdf"
#Hfile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Lung-n6_0.0_0.5_H.txt"
#outFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Lung-n6_0.0_0.5_H_plot.pdf"
#Hfile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Brain_Cerebellum-n6_0.0_0.5_H.txt"
#outFile <- "/scratch/users/paedugar/somaticMutationsProject/mutationSignatures/mutationSignaturesNMF/Brain_Cerebellum-n6_0.0_0.5_H_plot.pdf"
HEATMAP_GRADIENT <- colorRampPalette(c("#FFF5F0", "#3399ff"))(20)
HEATMAP_CELL_SIZE <- 0.1



#---------------------------
# METHODS

createHeatmap <- function(x, heatmapFile) {
    
    require(gplots)
    x <- as.matrix(x)
    
    #----
    # Gathering pars for heatmap2
    heatPars <- list(x = x,
                     col = HEATMAP_GRADIENT, 
                     dendrogram = "none", Rowv = FALSE,
                     cexCol = 0.001,
                     xlab = "Individuals", ylab = "Signature",
                     key.title = "Signature proportion",
                     symbreaks = F,
                     #scale = "column",
                     trace = "none", density = "none"
                     )
    widthPlot <- ncol(x) * HEATMAP_CELL_SIZE
    heightPlot <- widthPlot / 2
    
    
    #----
    
    pdf(heatmapFile, height = heightPlot, width = widthPlot)
    p <- do.call(heatmap.2, heatPars)
    dev.off()
    
    return(list(heatmap = p))
    
}



#---------------------------
# MAIN

#####
# Reading mutation data
# Reading table
Hmat <- read.table(Hfile, header = T, sep = "\t", stringsAsFactors = T)

# Calculate frequencies
Hmat <- t(t(Hmat) / (colSums(Hmat)))

p <- createHeatmap(Hmat, outFile)

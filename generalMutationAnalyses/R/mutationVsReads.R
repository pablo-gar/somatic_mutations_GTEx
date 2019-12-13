# Takes all the mutation counts for a given tissue/MAF combination
# and the read count table
# Plots a scatter plot of mapped reads vs mutation number

# Usage
# Rscript mutationVsReads.R outScatter.pdf readCountTable.txt count1.txt [count2.txt] [...] 
   
if(!dir.exists("../../R"))
   stop("Can't find the share R path, make sure to run script from directory where it lives")

source("../../R/plots.R")
source("../../R/ggthemes.R")

# Read count table should be produced with ../bin/getReadCounts

args <- commandArgs(T)
outScatter <- args[1]
readCountFile <- args[2]
args <- args[3:length(args)]

#workingDir <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/count/Liver/n6_0.0_0.5/"
#args <- file.path(workingDir, list.files(workingDir))

# Read read counts
readCountTable <- read.table(readCountFile, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

# Read mutation count tables
countTable <- data.frame(sample = "", mutations = rep(0, length(args)), stringsAsFactors = F)

for (i in 1:length(args)) {
    mutation <- read.table(args[i], sep = "\t", stringsAsFactors = F, header = F)
    countTable[i, "sample"] <- gsub(".txt", "", basename(args[i]))
    countTable[i, "mutations"] <- sum(mutation)
}

countTable$mappedReads <- readCountTable[countTable$sample, "n_reads"]
countTable <- countTable[!is.na(countTable$mappedReads),]

if(nrow(countTable) > 1) {
    p <- scatter(countTable, "mappedReads", "mutations", regression = T, xlab = "Uniquely mapped reads", ylab = "Mutations")
    p <- p + theme_noGrid()
} else {
    p <- ggplot(data.frame(x = 1, y = 1), aes(x = x, y = y)) + geom_point()
}

ggsave(outScatter, p, width = 5, height = 5)


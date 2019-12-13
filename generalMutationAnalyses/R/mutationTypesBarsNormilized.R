library("ggplot2")
library("reshape")
nucleotideFile <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/genome/Hg19_UCSC_knownGenes_exons_notOverlaping.fasta_oligoFreq_n1"
nucleotideCounts <- read.table(nucleotideFile, sep = "\t", stringsAsFactors = F, header = F)
rownames(nucleotideCounts) <- nucleotideCounts[,1]

barColors <- c("#3BC9F3", "#2B2E34", "#FC3218", "#CAD0CE", "#9CD169", "#F1CAC9")

args <- commandArgs(T)

outplot <- args[1]
args <- args[-1]

# Read tables
countTable <- list()
for (i in 1:length(args)) {
    mutation <- read.table(args[i], sep = "\t", stringsAsFactors = F, header = F)
    colnames(mutation) <- c("A", "T", "G", "C")
    mutation$from <- c("T", "C")
    
    mutation <- melt.data.frame(mutation)
    mutation$mut <- paste0(mutation[,1], ">", mutation[,2])
    
    # Normilize by nucleotide counts
    mutation$value <- mutation$value / (nucleotideCounts[mutation$from, 2] * 2)
    
    # Drop individual bases
    mutation <- mutation[,-c(1,2)]
    
    countTable[[i]] <- mutation
}
countTable <- do.call(rbind, countTable)

# Calculate average
statsMean <- tapply(countTable$value, countTable$mut, mean)
statsSD <- tapply(countTable$value, countTable$mut, sd)
stats <- data.frame(mut = names(statsMean), mean = statsMean, sd = statsSD, stringsAsFactors = F)
stats <- stats[!stats$mut %in% c("C>C", "T>T"), ]
stats$mut <- factor(stats$mut, levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), ordered = T)


p <- ggplot(stats, aes(x = mut, y = mean, fill = mut)) +
geom_bar(stat = "identity") + 
geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1 ) + 
scale_fill_manual(values = barColors) +
theme_bw()

ggsave(outplot, p)

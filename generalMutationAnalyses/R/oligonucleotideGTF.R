# Do module load mariadb/10.2.11

#' Gets the k-mer counts of all transcripts in a gtf files
#' @param size of k-mers
#' @param path to gtf file
#' @param path to genome fasta file (must be indexed)
#' @param path to output table, first column is transcript id, further columns are the counts 

library("GenomicFeatures")
library("Rsamtools")

args <- commandArgs(T)
oligoSize <- as.numeric(args[1]) #1
gtfFile <- args[2] #"/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/gencode.v19.genes.v7.patched_contigs.gtf"
fastaFile <- args[3] #"~/data/humanGenome/GRCh37.p13.maskedGenome.patched.fasta"
outFile <- args[4]

# Read gtf file coordinates
geneDB <- makeTxDbFromGFF(gtfFile, format = "gtf")
allTrans <- transcripts(geneDB)
allTrans <- renameSeqlevels(allTrans, paste0("chr", as.character(seqlevels(allTrans))))
allTrans <- allTrans[!grepl("MT", seqnames(allTrans))]

# Opens fasta file
fa <- open(FaFile(fastaFile))

allOligo <- list()
geneN <- 500
j <- 1
for(i in seq(1, length(allTrans), geneN)) {
    cat("Working with genes ", i, " - ", i + geneN - 1 , "\n")
    
    # Get gtf subset
    currentTrans <- allTrans[i:ifelse( i + geneN - 1 < length(allTrans), i + geneN - 1, length(allTrans))]
    # Extract sequence
    DNAstrings <- scanFa(fa, currentTrans)
    # Count oligos
    allOligo[[j]] <- oligonucleotideFrequency(DNAstrings, oligoSize)
    j <- j +1
}

allOligo <- as.data.frame(do.call(rbind, allOligo))
allOligo$gene_id <- allTrans$tx_name
allOligo <- allOligo[,c(ncol(allOligo), 1:(ncol(allOligo) - 1))]

# Saves file
write.table(allOligo, outFile, sep = "\t", quote = F, col.names = T, row.names = F)



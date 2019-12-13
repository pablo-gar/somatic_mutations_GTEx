# Given a fasta file gets the counts of all possible oligonucleotides
# of length n

require("Biostrings")


args <- commandArgs(T)
n <- as.numeric(args[1])
fastaFile <- args[2]
outFile <- args[3]

fastaSeq <- readDNAStringSet(fastaFile)
oligoFreq <- colSums(oligonucleotideFrequency(fastaSeq, n))


write.table(as.data.frame(oligoFreq), outFile, sep = "\t", col.names = F, row.names = T, quote = F)

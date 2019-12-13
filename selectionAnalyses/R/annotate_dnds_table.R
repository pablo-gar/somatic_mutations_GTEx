library(biomaRt)

args <- commandArgs(T)

mutationFile <- args[1] # sel_cv.txt from dndsout script
outputAnnotatedFile <- args[2] # 

#mutationFile <- "/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood/n6_0.0_0.5/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_sel_cv.txt"

# Getting ensembl ready
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Reading data
mutationTable <- read.table(mutationFile, sep = "\t", stringsAsFactors = F, header = T)
mutationTable <- mutationTable[!(mutationTable$n_syn == 0 & mutationTable$n_mis == 0),]
rownames(mutationTable) <- mutationTable$gene_name

geneDesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = mutationTable$gene_name, mart =ensembl)
geneDesc <- geneDesc[!duplicated(geneDesc$external_gene_name),]
rownames(geneDesc) <- geneDesc$external_gene_name

mutationTable$gene_description <- NA
mutationTable[rownames(geneDesc), "gene_description"] <- geneDesc$description

# Look for different categories of genes


# Genes with highest number of mutations
n <- 50
interestRows <- rownames(head(mutationTable[ order(-mutationTable$n_mis),], n))
mutationTable$highestMissSense <- 0
mutationTable[interestRows, "highestMissSense"] <- n:1

# Genes with lowest number of mutations
n <- 50
interestRows <- rownames(head(mutationTable[ order(mutationTable$n_mis, -mutationTable$n_syn),], n))
mutationTable$lowestMissSense <- 0
mutationTable[interestRows, "lowestMissSense"] <- n:1


# Genes with lowest dn/ds ratios
n <- 50
interestRows <- rownames(head(mutationTable[ order(mutationTable$wmis_cv, -mutationTable$n_syn),], n))
mutationTable$lowestWmis <- 0
mutationTable[interestRows, "lowestWmis"] <- n:1

# Genes with highest dn/ds ratios
n <- 50
interestRows <- rownames(head(mutationTable[ order(-mutationTable$wmis_cv, -mutationTable$n_mis),], n))
mutationTable$highestWmis <- 0
mutationTable[interestRows, "highestWmis"] <- n:1


write.table(mutationTable, outputAnnotatedFile, sep = "\t", quote = F, col.names = T, row.names = F)

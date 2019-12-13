## Finds several artifacts
source("../R/mutationGeneAnnotation.R", chdir = T)

# Finding artifacts in mutations accross genes
# Requires to run the selection snakemake

mutations <- read.table("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood_EXO/n6_0.0_0.5/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_genesAllMutations.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

globalInds <- read.table("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood_EXO/n6_0.0_0.5/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_globalIdnds.txt",
                         sep = "\t", header = T, stringsAsFactors = F)

sel_cv <- read.table("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood_EXO/n6_0.0_0.5/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_sel_cv.txt",
                     sep = "\t", header = T, stringsAsFactors = F)

sel_loc <- read.table("/scratch/users/paedugar/somaticMutationsProject/selection/dndsout/Whole_Blood_EXO/n6_0.0_0.5/dndsout_,geneExpression:AllExpressed,selectionMaf:0%_100%,_sel_loc.txt",
                     sep = "\t", header = T, stringsAsFactors = F)

#############################
# Finding the gene with highest number of synonymous mutations
syn <- mutations[ mutations$impact == "Synonymous",]
mis <- mutations[ mutations$impact == "Missense", ]

geneCountSyn <- as.data.frame(table(syn$gene))
colnames(geneCountSyn) <- c("gene_name", "Synonymous")

geneCountMis <- as.data.frame(table(mis$gene))
colnames(geneCountMis) <- c("gene_name", "Missense")

# Merging counts
geneCount <- merge(geneCountSyn, geneCountMis, all = T)
geneCount[is.na(geneCount)] <- 0

# Append cds length
geneCount$cds_length <- queryRefCDS(geneCount$gene_name, "gene_name", "CDS_length")

# Normalizing mutations
geneCount$SynonymousNorm <- geneCount$Synonymous / geneCount$cds_length
geneCount$MissenseNorm <- geneCount$Missense / geneCount$cds_length

geneCount <- geneCount[ order(geneCount$SynonymousNorm, decreasing =T), ]

# Looking at histone clusters for artifacts, it has >> Syn and << miss
mutArt1 <- mutations[mutations$gene == "HIST1H2AI",]
mutArt1 <- mutArt1[order(mutArt1$pos, mutArt1$gene), ]

# Look for positions that have > 1 mutation per sample
weirdMutations <- paste0(mutations$chr, "_", mutations$pos, "_", mutations$sampleID)
weirdMutations <- table(weirdMutations)
weirdMutations <- weirdMutations[weirdMutations>1]
weirdMutations <- sort(weirdMutations)

head(geneCount)


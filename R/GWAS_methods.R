binGenome <- function(genomeFile, binSize = 1e6, exclude = NULL) {
    # Creates a data frame with ranges of binned genome
    # genomeFile - string - path to a file with chromosome sizes, 
    #   tab delimitated, 2 columns, chromosome name and size on bps
    # bin.size - numeric - size of the bins

    # Reads file
    genome <- read.table(genomeFile, sep = "\t", stringsAsFactors = F, header = F)


    binResults <- list()

    for (i in 1:nrow(genome)) {
        chromName <- genome[i,1]
        chromName <- gsub("chr", "", chromName)
        if (sum(chromName %in% exclude) > 0) next
        chromSize <- genome[i,2]
        nWindows <- ceiling( chromSize / binSize)
        start <- 1
        end <- start + binSize
        for (j in 1:nWindows) {
            if (end > chromSize)
                end <- chromSize
            binResults <- c(binResults, list(c(chr = chromName, start = start, end = end)))
            start <- end + 1
            end <- start + binSize
        }
    }

    return (binResults)
}

lapplyHelper <- function(genomeRange, vcfFile, phenotype, genotype, topNsignificant, covariateMatrix = NULL, minorAlleleFreq, permute = T) {
    # This is a wrapper to use in a lapply function to call a GWAS call
    # on as single genomeRanges
    # genomeRange - list/data.frame dim(1,3) - has the following names {chr, start, end}
    
    results <- GWAS_Matrix(vcfFile = vcfFile, phenotype = phenotype, yieldSize = NULL,
                           chrom = genomeRange[1], startB = as.numeric(genomeRange[2]), endB = as.numeric(genomeRange[3]),
                           genotype = genotype, topNsignificant = topNsignificant, covariateMatrix = covariateMatrix, minorAlleleFreq = minorAlleleFreq,
                           permute = permute)
    return(results)
}

mergeGWASresults <- function(resultList, topNsignificant) {

    for(i in resultList){

        # Converting to chars the first to columns
        # I need to this or the merging below won't work
        i$realTop$snps <- as.character( i$realTop$snps )
        i$realTop$gene <- as.character( i$realTop$gene )

        # Creating final result object
        if(!exists("GWASresults")) {
            GWASresults <- i
            next
        }

        # Compiling list of top results
        # allocating memory
        temp1 <- i$realTop[ rep(1, 2*topNsignificant), ]

        # Appding the new with the previous 
        temp1[1:topNsignificant,] <- i$realTop
        temp1[(topNsignificant + 1) : (2*topNsignificant), ] <- GWASresults$realTop

        #And selecting the top hits
        temp1 <- temp1[order(temp1$pvalue),]
        GWASresults$realTop <- head(temp1, topNsignificant)


        # Compiling lowest pvalue from permutation
        GWASresults$permutedMin <- pmin(GWASresults$permutedMin, i$permutedMin)

        # Summing number of tests done
        GWASresults$totalSNPs <- GWASresults$totalSNPs + i$totalSNPs

    }

    return(GWASresults)


}



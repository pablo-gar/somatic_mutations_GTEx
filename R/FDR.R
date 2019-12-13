######################################
## Analyzes the results of GWASbyGene.r
#
# Example:
# Rscript resultStats.r ~/scripts/postSNVanalyses/GWASglobalSusceptibility/GWASbyGene/results.RData /scratch/users/paedugar/uvProject/SNVfiles/mutCount/geneMaps/phenoMatrix/phenoMatrix_relativeCountWithinRange_C\>T_FC.txt


#source("~/scripts/rModules/GWAS.r")
#source("~/scripts/rModules/fileManagement.r")
#source("~/scripts/rModules/superApply.r")

library("dplyr")

########
## METHODS


## CALCULATING FDRS AT DIFFERETN PVALUES FOR TWO SETS OF PVALUES, REAL AND PERMUTAED

getLowerTail <- function(x, values) {
	
	# Gets the sum of elements in x that are less than each element in values
	# Can handle repeated and unsorted elements in values
	
	# Stores original info for latter rearrangment
	original <- values
	
	# Sorts and unique
	values <- unique(values)
	sortedValues <- sort(values)
	sortedValues <- c(-Inf, sortedValues, Inf)
	
	# Makes the counting
	counts <-  cumsum(table(findInterval(c(sortedValues, x),sortedValues) ) - 1)
	# Ignonres last element
	counts <- counts[1:length(values)]
	
	# Reorder according to original vector
	names(counts) <- order(values)
	counts <- counts[sort(as.numeric(names((counts))))]
	names(counts) <- as.character(values)
	
	counts <- counts[as.character(original)]
	
	return(counts)
	
}

FDRcalculation <- function(realPvals, permPvals) {
	
	# Gets percentage of permPvals that less than each realPvalue
	expectedPercentage <- getLowerTail(permPvals, realPvals) / length(permPvals)
	
	# Gets expected number of positive (based on permutations) for each pvalue
	expectedFalsePos <- expectedPercentage * length(realPvals)
	
	# Calculate FDR, falsePos / real Pos
	FDR <- expectedFalsePos/(1:length(realPvals))

	# Re order according to original vector
	names(FDR) <- order(realPvals)
	FDR <- FDR [sort(as.numeric(names(FDR)))]
	names(FDR) <- as.character(realPvals)
	
	return(FDR)
}

#getFDRMatQTL <- function(matQTL, permDir, FDR_cutoff = 0.02) {
#	
#	# This a specific function only used in the GWASbyGene pipeline
#	# matQTL is a single reseult of matrixQTL containing ALL pvalues
#	# and an extra slot with the id of the gene
#	# permDir is the directory where permutations of each gene are saved
#	# as individual RData files this files were created with=
#	# lmMatGRangePhenomatrixDeletme from the GWAS pipeline
#	
#	gene <- matQTL$geneName
#	permFile <-joinPath(permDir, paste0(gene, ".RData")) 	
#	
#	if(file.exists(permFile)) {
#		load(permFile)
#	} else {
#		return()
#	}
#	
#	real <- matQTL$all$eqtls
#	perm <- lmresults$all$eqtls$pvalue
#	
#	#browser()
#	FDR <- FDRcalculation(real$pvalue, perm)
#	
#	positiveCount <- FDR <= FDR_cutoff
#	if(sum(positiveCount) > 0) {
#		real <- real[positiveCount,]
#		real$FDR <- FDR[positiveCount]
#		return(real)
#	} else {
#		return()
#	}
#	
#}
	
# FDR ESTIMATION BASED ON A SINGLE VECTOR OF PVALUES

getFDRvector <- function(pvals) {
	
	# Estimates FDR for each pvalue in the vector based on the expectation
	# of a uniform distribution of pvalues
	
	if(sum(pvals < 0 | pvals > 1) > 0)
		stop("Pvalues have to be between 0 and 1")
	
	# pvalues have to be sorted
	
	pvalUnique <- unique(sort(pvals))
	pvalCounts <- table(pvals)
	pvalCounts <- cumsum(pvalCounts)
	
	expectedItems <- length(pvals)*pvalUnique
	names(expectedItems) <- as.character(pvalUnique)
	
	#browser()
	FDRs <- expectedItems / pvalCounts
	FDRs[ FDRs > 1 ]  <- 1
	
	names(FDRs) <- as.character(names(expectedItems))
	FDRs <- FDRs[as.character(pvals)]
	names(FDRs) <- names(pvals)
	
	
	return(FDRs)
	
}

plotPvalFDR <- function(pvals, poinType = 19, pointSize = 0.3, ...) {
	
	FDR <- getFDRvector(pvals)
	plot (sort(pvals), FDR[ order(pvals) ], type = "l", xlab = "p-value", ylab = "FDR", ...)
	points(sort(pvals), FDR[ order(pvals) ], pch = poinType, cex = pointSize)
	
}


#####
# Other functions

#' 1. Calculates correlations between x and y column grouping by group_col
#' 2. Within group permutes y n_perm times and re-calculates correlations 
#' 3. Gets FDR for wilcox.comparisons between real and permutations at a p 
#'    (e.g if p = 0.05, then FDR based on a 5% chance of observing a significant test)
#' 
#' @param dataframe - a tidy object with x, y columns as wall a grouping column
#' @param x - column from dataframe in dplyr format 
#' @param y - column from dataframe in dplyr format 
#' @param group_col - column from dataframe in dplyr format 
#' @param p - p value cutoff for wilcox tests
#' @param n_perm - number of permutations
#' @param method - method to use for correlations, see ?cor.test 


cor_by_group_fdr <- function(dataframe, x, y, group_col, n_perm, p, method = "spearman", theoretical_fdr = T) { 
    
    x <- enquo(x)
    y <- enquo(y)
    group_col <- enquo(group_col)
    
    
    # Real pvalue
    real_cor <- dataframe %>%
        group_by(!!group_col) %>%
        summarise(cor_result = cor_helper(!!x, !!y, permute = F, method = method)) %>%
        #summarise(n = n()) %>%
        ungroup() 
    
    real_p <- rep(1, n_perm)
    for(i in seq_along(real_p)) {
        perm_cor <- dataframe %>%
            group_by(!!group_col) %>%
            summarise(cor_result = cor_helper(!!x, !!y, permute = T, method = method)) %>%
            ungroup() 
        
        test_result <- wilcox.test(real_cor$cor_result, perm_cor$cor_result)
        
        real_p[i] <- test_result[["p.value"]]
    }
    
    if(theoretical_fdr) {
        fdr <-(n_perm * p) / sum(real_p <= p)
    } else {
        
        perm_p <- rep(1, n_perm)
        
        for(i in seq_along(perm_p)) {
            # create two permuted vectors each time
            real_cor <- dataframe %>%
                group_by(!!group_col) %>%
                summarise(cor_result = cor_helper(!!x, !!y, permute = T, method = method)) %>%
                ungroup() 
            
            perm_cor <- dataframe %>%
                group_by(!!group_col) %>%
                summarise(cor_result = cor_helper(!!x, !!y, permute = T, method = method)) %>%
                ungroup() 
            
            test_result <- wilcox.test(real_cor$cor_result, perm_cor$cor_result)
            
            perm_p[i] <- test_result[["p.value"]]
        }
        fdr <-  sum(perm_p <= p) / sum(real_p <= p)
    }
    
    return(fdr)
}

cor_helper <- function(x,y, permute = F, method = "sp") {
    
    if(permute)
        y <- y[sample(length(y), replace = T)]
    
    cor.test(x,y, method = method)[["estimate"]]
}

source("~/scripts/Rmods/superApply.r", chdir = T)
source("~/scripts/Rmods/fileManagement.r", chdir = T)

source("~/scripts/FraserLab/skinGTExProject/Rmods/GTEX.r")

GWAS_Matrix <- function(vcfFile, phenotype, yieldSize = 10e5, genomicRange = NULL, chrom = NULL, startB = NULL, endB = NULL, genotype = "DS", topNsignificant = 50, covariateMatrix = NULL, granges = NULL, minorAlleleFreq = 0.1) {
    
    # Description: Performs GWAS genome-wide or on a specified region. If phenotype is permuted phenovector then this will retrieve
    #   the best pvalues which can be use to calculate and FDR
    #
    # Parameters:
    #   - vcfFile - string vector of length 1 - filepath to genotype vcf file
    #   - phenotype - numeric vector of matrix - contains the phenotype values with column names as infividual names corresponding to the vcf file 
    #   - yieldSize - numeric - number of SNPs to load in memory per iteration
    #   - chrom - string - desired chromosome if there is only one region of interest (ignored if startB and endB absent)
    #   - startB - numeric - start position  if there is only one region of interest (ignored if region and endB absent)
    #   - endB - numeric - start position  if there is only one region of interest (ignored if region and startB absent)
    #   - genotype - string - the type of genotype to be used (DS=dosage; GL=?)
    
    require(VariantAnnotation)
    
    # Gets individual names
    if (!is.null(dim(phenotype))) {
        individuals <- colnames(phenotype)
        phenotype <- as.matrix(phenotype)
    } else {
        individuals <- names(phenotype)
        phenotype <- t(as.matrix(phenotype))
    }
    individuals <- findNames(vcfFile, individuals)
    
    # Prepares parameters to read vcf file
    if(!is.null(chrom) & !is.null(startB) & !is.null(endB)){
        param <-  GRanges(chrom, IRanges(startB, endB)) 
        tab <- TabixFile(vcfFile, yieldSize=NA_integer_)
    } else if (!is.null(granges)){
        param <- genomicRange
        tab <- TabixFile(vcfFile, yieldSize=NA_integer_)
    }else{
        stop("You have to specify a Range")
    }
    
    
    # Select only individuals for which we have genotypes
    
    # Iteration on vcfFile
    open(tab)
    
    GWASresults <- list(realTop = data.frame(snps = rep("", topNsignificant) , gene = "", statistic = 0, pvalue = 1, FDR = 1, beta = 0),
                        permutedMin = rep(1,nrow(phenotype)), 
                        totalSNPs = 0
                        )
    
    #count = 0
    while (nrow(vcf_yield <- readGeno(tab, genotype, param=param))){ 
        
        doRegression <- T
        
        #if (count > 1 ) break
        #count = count + 1
        
        firstPos <- rownames(vcf_yield)[1]
        lastPos <- rownames(vcf_yield)[nrow(vcf_yield)]
    
        cat("Analyzing from: ", firstPos, "   to: ", lastPos, "\n")
        
        # Get best imputed genotype
        cat("\t", as.character(Sys.time()), " Getting genotypes\n")         
        
        if(genotype == "DS"){
            genoMat <- vcf_yield            
        }else if (genotype == "GT") {
            genoMat <- getGenoFromGT(vcf_yield)
        }else if (genotype == "GL") {
            genoMat <- getGeno(vcf_yield)
        } else {
            stop ("Unrecognized genotype mode")
        }
        
        
        # Checking whether we have correct genotypes
        if(is.null(nrow(genoMat))) 
            doRegression <- F
        if(doRegression & nrow(genoMat) < 1) 
            doRegression <- F 
        if(is.null (colnames(genoMat)))
            colnames(genoMat) <- individuals
            
        
        
        # Selecting individuals present in both the genotype and phenotype vector
        # Selecting also SNPs that pass the MAF filter
        if(doRegression) { 
            
            # Arraging data and selecting only individuals for which we have phenotypes
            genoMat <- as.matrix(genoMat)  
            colnames(genoMat) <- gsub("(GTEX-.+?)-.+", "\\1", colnames(genoMat))
            
            phenotype <- phenotype[,colnames(genoMat)]
            if (!is.null(covariateMatrix)) 
                covariateMatrix <- covariateMatrix[,colnames(genoMat)]
            
            
            # Only positions that have the MAF requirement
            cat("\t", as.character(Sys.time()), " Getting positions with variable genotypes (MAF = ", minorAlleleFreq, ")\n")     
            indeces <- apply(genoMat, 1, selectMAF, minorAlleleFreq = minorAlleleFreq)
            
            # Only covariates that have at least 3 variable individuals
            indeces <- apply(covariateMatrix, 1, function(x) { min(table(x)) < length(x) - 3})
            covariateMatrix <- covariateMatrix[indeces, , drop = F]
            
            # Don't do regression if MAF threshold is not met
            if(sum(indeces) == 0) doRegression <- F
        }
        
        # If all checks pass do regression
        if(doRegression) { 
        
            genoMat <- genoMat[indeces, , drop = F]

            cat("\t", as.character(Sys.time()), " Calling linear regression against phenotype\n")                   
            
                #Recording real results data
                realResults <- lmMatFunction (genoMat, phenotype[1, , drop = F], useModel = modelLINEAR, coovariates = character(), pvalCutoff = max(GWASresults$realTop$pvalue), cvrt = covariateMatrix,  min.pv.by.genesnp = F, noFDRsaveMemory = F)
                realResults <- rbind(realResults$all$eqtls, GWASresults$realTop)
                realResults <- realResults[order(realResults$pvalue), ]
            
                #Do permutations (keep only the minumun so far for each permutation)
                permutedResults <- lmMatFunction (genoMat, phenotype, useModel = modelLINEAR, coovariates = character(), pvalCutoff = 1, cvrt = covariateMatrix, min.pv.by.genesnp = T, noFDRsaveMemory = TRUE)
                
                GWASresults$realTop <- head(realResults, n = topNsignificant)
                GWASresults$permutedMin <- pmin(GWASresults$permutedMin, permutedResults$all$min.pv.gene)
                GWASresults$totalSNPs = GWASresults$totalSNPs + nrow(genoMat)
                #GWASresults$geno = genoMat
                #GWASresults$pheno = phenotype
                
                
            
            gc();gc()
        }
        
        cat("\t", as.character(Sys.time()), " Done\n")                  
        
        #Exit if not yieldsize
        if(!is.null(chrom) & !is.null(startB) & !is.null(endB))
            break
        
        #Exit if granges as input
        if(!is.null(granges))
            break

        
    }
    close(tab)
    return(GWASresults)
}
    
    



GWAS2 <- function(vcfFile, phenotype, yieldSize = 10e5, tasks, workingDir, chrom = NULL, startB = NULL, endB = NULL, dosage = T,extraScriptLines){
    
    require(VariantAnnotation)
    
    # Reads phenotypes and individuals
    phenoInd <- readPhenotype(vcfFile, phenotype)
    phenotype <- phenoInd[[1]]
    individuals <- phenoInd[[2]]
    

    # Making iterator
    if (dosage){
        genoType <- "DS"
    }else{
        genoType <- "GL"
    }
    
    if(!is.null(chrom) & !is.null(startB) & !is.null(endB)){
        param <- ScanVcfParam(fixed="ALT", geno=genoType, samples = individuals, which = GRanges(chrom, IRanges(startB, endB)) )
        tab <- TabixFile(vcfFile, yieldSize=NA_integer_)
    }else{
        if(is.na(yieldSize))
            stop("yieldSize not found")
        if(yieldSize < 2)
            stop("yieldsize has to be greated than 2")
        
        param <- ScanVcfParam(fixed="ALT", geno=genoType, samples = individuals)    
        tab <- TabixFile(vcfFile, yieldSize=yieldSize)
    }
    
    
    # Iterating over vcfFile
    GWASresults <- list()
    open(tab)
    while (nrow(vcf_yield <- readGeno(tab, genoType, param=param))){
        firstPos <- rownames(vcf_yield)[1]
        lastPos <- rownames(vcf_yield)[nrow(vcf_yield)]
    
        cat("Analyzing from: ", firstPos, "   to: ", lastPos, "\n")
        
        # Get best imputed genotype
        cat("\t", as.character(Sys.time()), " Getting genotypes\n")         
        
        if(dosage){
            genoMat <- vcf_yield            
        }else{
            genoMat <- getGeno(vcf_yield)
        }
        
        genoMat <- as.data.frame(t(genoMat)) # Positions as columns, inds as rows
        
        
        # Only positions that have at least one ind with a different phenotype
        cat("\t", as.character(Sys.time()), " Getting positions with variable genotypes\n")         
        indeces <- sapply(genoMat, function(x) length(unique(x)) > 1)
        genoMat <- genoMat[,indeces]

        cat("\t", as.character(Sys.time()), " Calling linear regression against phenotype\n")                   
            tempR <- superApply(genoMat, linearRegressionHelper,  y = phenotype, tasks = tasks, workingDir = workingDir, clean = F, extraScriptLines = extraScriptLines)
        GWASresults <- c(GWASresults, tempR)
        
        gc();gc()
        cat("\t", as.character(Sys.time()), " Done\n")                  
        
        #Exit if not yieldsize
        if(!is.null(chrom) & !is.null(startB) & !is.null(endB)){
            break
        }
        
    }
    close(tab)
    
    return(do.call(c,GWASresults))
    
}


# Performs a GWAS on a vectors whose names are GTEX-individuals
# Takes a vcfFile as input

GWAS <- function(vcfFile, phenotype, yieldSize = 10e3, parallel = F, tasks = NULL, chrom = NULL, startB = NULL, endB = NULL, dosage = F){
    
    require(VariantAnnotation)
    
    if(parallel){
        require(BiocParallel)
        if(is.null(tasks))
            stop("when parallel is true, task has to be specified")

        mParam <- MulticoreParam(workers = tasks)
    } 
    
    # Reads phenotypes and individuals
    phenoInd <- readPhenotype(vcfFile, phenotype)
    phenotype <- phenoInd[[1]]
    individuals <- phenoInd[[2]]
    

    # Making iterator
    if (dosage){
        genoType <- "DS"
    }else{
        genoType <- "GL"
    }
    
    if(!is.null(chrom) & !is.null(startB) & !is.null(endB)){
        param <- ScanVcfParam(fixed="ALT", geno=genoType, samples = individuals, which = GRanges(chrom, IRanges(startB, endB)) )
    }else{
        if(is.na(yieldSize))
            stop("yieldSize not found")
        if(yieldSize < 2)
            stop("yieldsize has to be greated than 2")
        
        param <- ScanVcfParam(fixed="ALT", geno=genoType, samples = individuals)    
    }
    
    tab <- TabixFile(vcfFile, yieldSize=yieldSize)
    
    # Iterating over vcfFile
    GWASresults <- list()
    open(tab)
    while (nrow(vcf_yield <- readGeno(tab, genoType, param=param))){
        firstPos <- rownames(vcf_yield)[1]
        lastPos <- rownames(vcf_yield)[nrow(vcf_yield)]
    
        cat("Analyzing from: ", firstPos, "   to: ", lastPos, "\n")
        
        # Get best imputed genotype
        cat("\t", as.character(Sys.time()), " Getting genotypes\n")         
        
        if(dosage){
            genoMat <- vcf_yield            
        }else{
            genoMat <- getGeno(vcf_yield)
        }
        
        genoMat <- as.data.frame(t(genoMat)) # Positions as columns, inds as rows
        
        
        # Only positions that have at least one ind with a different phenotype
        cat("\t", as.character(Sys.time()), " Getting positions with variable genotypes\n")         
        indeces <- sapply(genoMat, function(x) length(unique(x)) > 1)
        genoMat <- genoMat[,indeces]

        cat("\t", as.character(Sys.time()), " Calling linear regression against phenotype\n")                   
        if(parallel){
            tempR <- bplapply(genoMat, linearRegressionHelper, y = phenotype, BPPARAM = mParam)
        }else{
            tempR <- lapply(genoMat, linearRegressionHelper, y = phenotype)
        }
        GWASresults <- c(GWASresults, tempR)
        
        gc();gc()
        cat("\t", as.character(Sys.time()), " Done\n")                  
        
        #Exit if not yieldsize
        if(!is.null(chrom) & !is.null(startB) & !is.null(endB)){
            break
        }
        
    }
    close(tab)
    
    return(do.call(c,GWASresults))
    
}

#######
## FUNCTIONS FOR GWAS BY GENE

## Gets all SNPs in a given region from a VCF file
regionSNPs <- function(vcfFile, individuals, chrom, startB, endB, grange = NULL, genotype = "DS", chromFormat = "ENCODE"){ # Do not modify genotype or else function will not workd
    require(VariantAnnotation)
    require(GenomicRanges)
    #Transforming chrom name format
    if(!chromFormat %in% c("ENCODE", "ENSEMBL")) stop ("Invalid chromosome format, options: [ENCODE|ENSEMBL]")
    
    # R Reading paramters and genotype table
    individuals <- findNames(vcfFile, individuals)
    
    # Either explicit coordinates or granges
    scanParams <- list(fixed="ALT", geno=genotype, samples = individuals)
    if(is.null(grange)){
        if (chromFormat == "ENCODE") chrom <- gsub("chr", "", chrom)
        scanParams<- c(scanParams, list(which = GRanges(chrom, IRanges(startB, endB))) )
    }else{
        if (chromFormat == "ENCODE") {seqlevels(grange) <- gsub("chr", "", seqlevels(grange));seqnames(grange) <- gsub("chr", "", seqnames(grange))}
        scanParams<-  c(scanParams,list(which = grange)) 
    }
    
    
    param <- do.call(ScanVcfParam,  scanParams)
    genoMat <- readGeno(vcfFile, genotype, param=param)
    
    
    if(genotype == "GL")
        genoMat <- getGeno(genoMat)
    if(genotype == "GT")
        genoMat <- getGenoFromGT(genoMat)
    
    genoMat <- as.data.frame(t(genoMat))
    
    if (ncol(genoMat) < 1){
        warning("No SNPs found for range")  
        return(NULL)
    }
    
    # Need to add because for some strange reason when there is only one GRange it does not get the names
    if(length(scanParams$which) > 1){
        
        tempParams <- scanParams
        tempParams$which <- tempParams$which[1]
        
        
        param <- do.call(ScanVcfParam,  tempParams)
        individualNamesInGeno <- colnames(readGeno(vcfFile, genotype, param=param))
        rownames(genoMat) <- individualNamesInGeno
    }
        
    
    # Getting only varaible genotypes
    indeces <- sapply(genoMat, function(x) length(unique(x)) > 1)
    
    if (sum(indeces) < 1){
        warning("No SNPs with variable genotypes were found")   
        return(NULL)
    }
    
    # Saving names when only one SNP has been found
    flagNames = FALSE
    if(ncol(genoMat) == 1) {
        genoMatNames <- rownames(genoMat)
        flagNames <- TRUE
    }
    
    genoMat <- genoMat[,indeces, drop = F]
    #if(flagNames)
        #colnames(genoMat) <- genoMatNames
    
    return(genoMat)

}


##
# Given SNP ids with the GTEX format gets the genotypes of the given individuals from the given vcf file
# ids are of the form: 1_20301_[...]
getSNPsByID <- function(ids, vcfFile, individuals, genotype = "DS", chromFormat = "ENSEMBL") {
    ids <- strsplit(ids, "_")
    chrs <- sapply(ids, function(x) x[1])
    starts  <- as.integer(sapply(ids, function(x) x[2]))
    
    
    snps <- regionSNPs(vcfFile = vcfFile, individuals = individuals, chrom = chrs, startB = starts, endB = starts + 1, genotype = genotype, chromFormat = chromFormat)
    
    return(snps)
        
}

readSNPs <- function(ids, vcfFile, individuals, genotype, chromFormat) {
        
    # Reads in the snps genotypes of the given individuals and presents them in a
    # table format of the type individuals as columns and phenotypes as rows

    snpGeno <- getSNPsByID (ids = ids, vcfFile = vcfFile, individuals = individuals, genotype = genotype, chromFormat = chromFormat)
    rownames(snpGeno) <- GTEX$simplifyId(rownames(snpGeno))
    snpGeno <- snpGeno[individuals, , drop = F]
    snpGeno <- t(snpGeno)

    return(snpGeno)
    
}


# Looks for all SNPS in, upstream and downstream a single GRange
# if multiple GRagnes given only the first will be considered

rangeSNPs <- function(grange, upstream = 0, downstream = 0, vcfFile, individuals, genotype = "DS"){
    require(VariantAnnotation)
    require(GenomicRanges)
    # Expand GRagnes
    grange <- grange[1]
    grange <- expandGRanges (grange, upstream, downstream)
    # Getting SNPs
    return(regionSNPs(vcfFile, individuals, grange = grange, genotype = genotype))
}

# Helper function that expands granges upstream and downstream by a given size
expandGRanges <- function(grange, upstream = 0, downstream = 0){
    
    require(GenomicRanges)
    
    # Getting grange
    shiftUp <-ifelse(strand(grange) == "-", upstream, -upstream)
    
    #Resize GRange to include upstream and downstream regions
    grange <- resize(shift(grange, shiftUp), width(grange) + upstream + downstream)
    
    return(grange)
}


# Returns the snp location of given SNP ids
getSNPlocation <- function(ids) {
    
    require(biomaRt)
    
    SNP.DB <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="grch37.ensembl.org", path="/biomart/martservice")
    
    results <- getBM(c("refsnp_id","allele","chr_name","chrom_start"),  #, "chrom_strand","associated_gene", "ensembl_gene_stable_id"), 
                     filters="snp_filter", 
                     values=ids, 
                     mart=SNP.DB)
    return(results)
}

# Takes a gene name, a GRanges object containing the coordinates of the gene (whose names are geneNames),
# and a phenoMatrix (constructed with joinMaps_andCreateMutMatrixCluster.py) with rows as Gene names and 
# columns as individuals, and performs a linear regression for all SNPs in the specified region against the
# values in the phenoMatrix

lmFunction <- function(x,y, na.rm = F, inf.rm = F, coeff = F) {
    lapply(as.data.frame(x), linearRegressionHelper, y = y, na.rm = na.rm, inf.rm = inf.rm, coeff = coeff)
}

lmGRangePhenomatrix <- function(x, vcfFile, granges, phenoMatrix, upstream = 0, downstream = 0, permutations = NULL, snpIds = NULL){
    require(GenomicRanges)
    require(VariantAnnotation)
    
    #flush.console()
    #cat(x, "\n")
    
    if(! is.vector(x)) stop("x is not a vector")
    if( length(x) > 1 ) stop("Only one gene can be specified")
    if(!is.null(permutations) && !is.numeric(permutations)) stop("permutations has to be numeric and indicate the number of permutations desired")
    
    inds <- colnames(phenoMatrix)
    phenoVector <- phenoMatrix[x,]
    currentGene <- granges[x]
    
    #Getting SNPs
    if(is.null(snpIds)) {
        snps <- rangeSNPs(currentGene, upstream, downstream, vcfFile, inds)
    } else {
        snps <-  getSNPsByID(snpIds, vcfFile = vcfFiles, individuals = inds, genotype = "DS", chromFormat = "ENSEMBL") 
    }
    
    if(is.null(snps)) return(NA)
    if( is.null(nrow(snps)) ) return(NA)
    
    rownames(snps) <- gsub("(GTEX-.+?)-.+", "\\1", rownames(snps))
    
    #Getting phenotypes from individuals 
    phenoVector<- as.numeric(phenoMatrix[x, rownames(snps)])
    #return(phenoVector)
    #Performing linear regression 
    if(is.null(permutations)){
        lmresults <-do.call(c, lmFunction (x = snps, y = phenoVector, na.rm = T, inf.rm = T))
    }else{ # This will permute phenotypes 1 mi times and then count how many of those the most
        # significant assiociation exceeds the most significant association of the original data
        
        #real data
        lmresults <-do.call(c, lmFunction (x = snps, y = phenoVector, na.rm = T, inf.rm = T))
        
        #permutes
        lmresultsPerms <- lapply(rep(T, permutations), 
                    function(x, snps, phenoVector){ phenoVector <- sample(phenoVector)
                                    return(do.call(c, lmFunction(x = snps, y = phenoVector,  na.rm = T, inf.rm = T)))
                                }, snps = snps, phenoVector = phenoVector
                )
        
        #In how many permutations do we exceed pval?
        lmresults <- sapply (lmresultsPerms, function(x, real)  min(x, na.rm = T) < real, real = min(lmresults, na.rm = T))
        lmresults <- sum(lmresults)/permutations
    }
    gc()
    return (lmresults)
}

#Same as the one above but using matrix-based calculations (way faster) 

lmMatFunction <- function(x,y, useModel = modelLINEAR, coovariates = character(), pvalCutoff = 1, errorCovariance = numeric(), outFile = tempfile(), min.pv.by.genesnp = F, noFDRsaveMemory = F, cvrt = NULL){
    require(MatrixEQTL)
    
    xMat <- SlicedData$new()
    yMat <- SlicedData$new()
    cvrtMat <- SlicedData$new()
    
    xMat$CreateFromMatrix(x)
    yMat$CreateFromMatrix(y)
    
    if(!is.null(cvrt)) 
        cvrtMat$CreateFromMatrix(cvrt)
    
    
    results <- Matrix_eQTL_engine(snps = xMat, gene = yMat, cvrt = cvrtMat, output_file_name = outFile,
                      pvOutputThreshold = pvalCutoff, useModel = useModel, min.pv.by.genesnp = min.pv.by.genesnp, noFDRsaveMemory = noFDRsaveMemory)
    
    return(results)
}

# IT is used for apply, it takes the eQTLtable that focuses in specific snps
lmMatGRangePhenomatrix.eQTLtable <- function(x, vcfFile, granges, phenoMatrix, upstream, downstream, permutations, eQTLtable) {
    
    snpIds <- eQTLtable$variant_id[ grepl(x, eQTLtable$gene_id)]
    if (length(snpIds) < 1 ) return (NULL)
    results <- lmMatGRangePhenomatrix(x,  vcfFile=vcfFile, granges=gtexGenes, phenoMatrix=phenoM, upstream = ups, downstream = dos, permutations = permutations, snpIds = snpIds)
    
    return(results) 

}


lmMatGRangePhenomatrix <- function(x, vcfFile, granges, phenoMatrix, upstream = 0, downstream = 0, pvalCutoff = 0.05, bonferroni = T,  permutations = NULL, snpIds = NULL, cvrt = NULL){
    require(MatrixEQTL)
    require(GenomicRanges)
    require(VariantAnnotation)
    
    flush.console()
    cat(x, "\n")
    
    if(! is.vector(x)) stop("x is not a vector")
    if( length(x) > 1 ) stop("Only one gene can be specified")
    if(!is.null(permutations) && !is.numeric(permutations)) stop("permutations has to be numeric and indicate the number of permutations desired")
    
    currentGene <- granges[x]
    inds <- colnames(phenoMatrix)
    phenoVector <- do.call(c,phenoMatrix[x,])
    
    realValues <- !is.na(phenoVector) & !is.infinite(phenoVector)
    phenoVector <- phenoVector[realValues]
    inds <- inds[realValues]
    
    #Getting SNPs
    
    #Getting SNPs
    if(is.null(snpIds)) {
        snps <- rangeSNPs(currentGene, upstream, downstream, vcfFile, inds)
    } else {
        snps <-  getSNPsByID(snpIds, vcfFile = vcfFile, individuals = inds, genotype = "DS", chromFormat = "ENSEMBL") 
    }   
    
    if( is.null(nrow(snps)) ) return(NA)
    
    rownames(snps) <- gsub("(GTEX-.+?)-.+", "\\1", rownames(snps))
    indsWithGeno <- rownames(snps)
    
    #Getting phenotypes from individuals 
    phenoVector<- as.numeric(phenoMatrix[x, indsWithGeno])
    
    #Selecting individuals for covariates
    if(!is.null(cvrt))
        cvrt <- cvrt[indsWithGeno,]
    
    # Transforming to matrices for matrixeQTL
    snps <- t(as.matrix(snps))
    phenoVector <- t(as.matrix(phenoVector))
    if(!is.null(cvrt))
        cvrt <- t(as.matrix(cvrt))
    
    # Getting pval cutoff
    if (bonferroni) pvalCutoff <- pvalCutoff / nrow(snps)
    #cat(pvalCutoff)
    
    if(is.null(permutations)){
        lmresults <- lmMatFunction (x = snps , y = phenoVector, useModel = modelLINEAR, coovariates = character(), 
                        pvalCutoff = pvalCutoff, errorCovariance = numeric(), outFile = tempfile(), cvrt = cvrt)

    }else{ # This will permute phenotypes 1 mi times and then count how many of those the most
        # significant assiociation exceeds the most significant association of the original data
        
        #real data
        lmresults <- lmMatFunction (x = snps , y = phenoVector, useModel = modelLINEAR, coovariates = character(), 
                        pvalCutoff = pvalCutoff, errorCovariance = numeric(), outFile = tempfile(),
                        min.pv.by.genesnp = T, cvrt = cvrt)

        #Creates a matrix of phenovectors  where all of them are permuted in order
        phenoPermuted <- as.matrix(do.call(rbind, lapply(rep(T, permutations), function(x,y) sample(y), y = phenoVector[1,])))
        colnames(phenoPermuted) <- colnames(phenoVector)
        
        #Perfoms linear regressions
        lmresultsPerms <- lmMatFunction (x = snps , y = phenoPermuted, useModel = modelLINEAR, coovariates = character(), 
                        pvalCutoff = pvalCutoff, errorCovariance = numeric(), outFile = tempfile(),
                        min.pv.by.genesnp = T, cvrt = cvrt)
        
        #In how many permutations do we exceed pval
        lmresults <- sum(lmresultsPerms$all$min.pv.gene <= lmresults$all$min.pv.gene[1])
    }
    gc()
    return (lmresults)
}

# Gets the genotype of n strongest associations for a given gene and a phenoveoctor

getStrongestSNPs <- function(x, nStrongest = 1, vcfFile, granges, phenoMatrix, upstream = 0, downstream = 0, pvalCutoff = 0.05, bonferroni = F, snpIds = NULL){
    
    # Getting associations
    lmresults <- lmMatGRangePhenomatrix(x, vcfFile = vcfFile, granges = granges, phenoMatrix = phenoMatrix, upstream = upstream, downstream = downstream, pvalCutoff = pvalCutoff, bonferroni = bonferroni, snpIds = snpIds)
    
    # Getting top hit snps
    if( lmresults$all$neqtls < nStrongest )
        nStrongest <- lmresults$all$neqtls
    
    snpIds <- as.character(lmresults$all$eqtl$snps[1:nStrongest])
    if(length(snpIds) < 1 ) return (NULL)
    
    # Getting genotypes of snps
    inds <- colnames(phenoMatrix)
    snps <- getSNPsByID(snpIds, vcfFile = vcfFile, individuals = inds, genotype = "DS", chromFormat = "ENSEMBL")
    
    if(is.null(dim(snps))) {
        snps <- as.matrix(snps)
        colnames(snps) <- snpIds
    }
    snps <- snps[ , colnames(snps) %in% snpIds, drop = F]
    
    rownames(snps) <- gsub("(GTEX-.+?)-.+", "\\1", rownames(snps))
    
    #Getting phenotypes from individuals 
    phenoVector<- as.numeric(phenoMatrix[x, rownames(snps)])
    
    snpsPheno <- cbind(snps, phenoVector)
    
    return (list(pvalues = lmresults$all$eqtl[1:nStrongest, ], snpsPheno = snpsPheno))
                       
    
}



lmMatGRangePhenomatrixDeletme <- function(x, vcfFile, granges, phenoMatrix, upstream = 0, downstream = 0, pvalCutoff = 0.05, bonferroni = T,  permutations = NULL, saveToDir = NULL){
    require(MatrixEQTL)
    require(GenomicRanges)
    require(VariantAnnotation)
    
    flush.console()
    cat(x, "\n")
    
    if(! is.vector(x)) stop("x is not a vector")
    if( length(x) > 1 ) stop("Only one gene can be specified")
    if(!is.null(permutations) && !is.numeric(permutations)) stop("permutations has to be numeric and indicate the number of permutations desired")
    
    currentGene <- granges[x]
    inds <- colnames(phenoMatrix)
    phenoVector <- do.call(c,phenoMatrix[x,])
    
    realValues <- !is.na(phenoVector) & !is.infinite(phenoVector)
    phenoVector <- phenoVector[realValues]
    inds <- inds[realValues]
    
    #Getting SNPs
    snps <- rangeSNPs(currentGene, upstream, downstream, vcfFile, inds)
    if( is.null(nrow(snps)) ) return(NA)
    
    rownames(snps) <- gsub("(GTEX-.+?)-.+", "\\1", rownames(snps))
    
    #Getting phenotypes from individuals 
    phenoVector<- as.numeric(phenoMatrix[x, rownames(snps)])
    #return(phenoVector)
    
    # Transforming to matrices for matrixeQTL
    snps <- t(as.matrix(snps))
    phenoVector <- t(as.matrix(phenoVector))
    
    # Getting pval cutoff
    if (bonferroni) pvalCutoff <- pvalCutoff / nrow(snps)
    #cat(pvalCutoff)
    
    if(is.null(permutations)){
        lmresults <- lmMatFunction (x = snps , y = phenoVector, useModel = modelLINEAR, coovariates = character(), 
                        pvalCutoff = pvalCutoff, errorCovariance = numeric(), outFile = tempfile())
    } else {
        phenoMPermuted <- do.call(rbind, lapply(rep(T, permutations), function(x,y) y[sample(1:length(y))], y = phenoVector))
        lmresults <- lmMatFunction (x = snps , y = phenoMPermuted, useModel = modelLINEAR, coovariates = character(), 
                        pvalCutoff = pvalCutoff, errorCovariance = numeric(), outFile = tempfile())
    }
        

    gc()
    
    if (is.null(saveToDir)) {
        return (lmresults)
    } else {
        mkdirRecursive(saveToDir)
        save(lmresults, file = joinPath(saveToDir, paste0(x,".RData")))
        return (TRUE)
    }
    
}

lmMatGRangePhenomatrixDeletme2 <- function(x, vcfFile, granges, phenoMatrix, upstream = 0, downstream = 0, pvalCutoff = 0.05, bonferroni = T,  permutations = NULL){
    require(MatrixEQTL)
    require(GenomicRanges)
    require(VariantAnnotation)
    
    flush.console()
    cat(x, "\n")
    
    if(! is.vector(x)) stop("x is not a vector")
    if( length(x) > 1 ) stop("Only one gene can be specified")
    if(!is.null(permutations) && !is.numeric(permutations)) stop("permutations has to be numeric and indicate the number of permutations desired")
    
    currentGene <- granges[x]
    inds <- colnames(phenoMatrix)
    #phenoVector <- do.call(c,phenoMatrix[x,])
    #
    #realValues <- !is.na(phenoVector) & !is.infinite(phenoVector)
    #phenoVector <- phenoVector[realValues]
    #inds <- inds[realValues]
    
    #Getting SNPs
    snps <- rangeSNPs(currentGene, upstream, downstream, vcfFile, inds)
    if( is.null(nrow(snps)) ) return(NA)
    
    rownames(snps) <- gsub("(GTEX-.+?)-.+", "\\1", rownames(snps))
    
    #browser()
    #Getting phenotypes from individuals 
    phenoVector<- (phenoMatrix[, rownames(snps)])
    #return(phenoVector)
    
    # Transforming to matrices for matrixeQTL
    snps <- t(as.matrix(snps))
    phenoVector <- (as.matrix(phenoVector))
    
    
    # Getting pval cutoff
    if (bonferroni) pvalCutoff <- pvalCutoff / nrow(snps)
    #cat(pvalCutoff)
    
    if(is.null(permutations)){
        lmresults <- lmMatFunction (x = snps , y = phenoVector, useModel = modelLINEAR, coovariates = character(), 
                        pvalCutoff = pvalCutoff, errorCovariance = numeric(), outFile = tempfile())
    } else {
        phenoMPermuted <- do.call(rbind, lapply(rep(T, permutations), function(x,y) y[sample(1:length(y))], y = phenoVector))
        lmresults <- lmMatFunction (x = snps , y = phenoMPermuted, useModel = modelLINEAR, coovariates = character(), 
                        pvalCutoff = pvalCutoff, errorCovariance = numeric(), outFile = tempfile())
    }
        

    gc()
    return (lmresults)
}

## Linear regression permutation of single position
lmSNPpermutation <- function(vcfFile, phenotype, genotype, nIter = 10000, parallel = F, tasks, param = NULL){
    
    require(VariantAnnotation)
    if(is.null(param)){
        if(parallel){
            require(BiocParallel)
            if(is.null(tasks))
                stop("when parallel is true, task has to be specified")

            mParam <- MulticoreParam(workers = tasks)
        }
    }else{
        require(BiocParallel)
        require(Rmpi)
        #register(param)
        #mParam <- param
    }
    
    
    #Reads phenotypes and individuals
    phenoInd <- readPhenotype(vcfFile, phenotype)
    phenotype <- phenoInd[[1]]
    individuals <- phenoInd[[2]]
    
    # Reads genotype
    genotype <- corGenoPheno(vcfFile, phenotype, genotype)
    
    # Performs linear regression in original data
    originalP <- linearRegressionHelper(phenotype, genotype)
    
    # Perfoms lienar regression with randomized data
    if(is.null(param)){ 
        if(parallel){
            newP <- bplapply(1:nIter, linearRegressionHelperPerm, phenotype = phenotype, genotype = genotype, BPPARAM = mParam)
        }else{
            newP <- lapply(1:nIter, linearRegressionHelperPerm, phenotype = phenotype, genotype = genotype) 
        }
        return(c(originalP, do.call(c,newP)))
    }else{
        #newP <- bplapply(1:nIter, linearRegressionHelperPerm, phenotype = phenotype, genotype = genotype)
        return(list(phenotype, genotype))
    }   
    
}


# Given a position name (genotype) gets the actual genotypes
# for the individuals contained in the phenotype vector(whose names are individuals)
corGenoPheno <- function(vcfFile, phenotype, genotype, dosage = T){
    require(VariantAnnotation)
    
    # Reads phenotypes and individuals
    phenoInd <- readPhenotype(vcfFile, phenotype)
    phenotype <- phenoInd[[1]]
    individuals <- phenoInd[[2]]
    
    # Reads genotype 
    chrom <- gsub("^(\\w+?)_.+", "\\1", genotype )
    startGeno <- as.integer(gsub("^\\w+?_(\\w+?)_.+", "\\1", genotype))
    endGeno <- startGeno + 1
     
    # Making iterator
    if (dosage){
        genoType <- "DS"
    }else{
        genoType <- "GL"
    }

    param <- ScanVcfParam(fixed="ALT", geno=genoType, samples = individuals, which = GRanges(chrom, IRanges(startGeno, endGeno)) )      
    
    # Reading VCF
    vcfGeno <- readGeno(vcfFile, genoType, param=param)[1,]
    
    return(vcfGeno)
}

# Reads a phenotype table of that has the following columns, mutation Type, gtexIndividual, age, race and whatever data wants to be read
readPhenotypeTable <- function(phenotypeTable, mut = "C>T",  gtexId= 1, mutId = 2, dataId = 3, age = 40 , race = 3, ageMax = 200) {

    phenoTable <- read.table(phenotypeTable, sep = "\t", stringsAsFactors = F, header = T)


    if(!mut %in% phenoTable[,mutId])
        stop("Mutation type not found in phenotype table")
    
    phenoTable <- phenoTable[ phenoTable[,mutId] == mut, ]
    phenoTable <- phenoTable[phenoTable[,"age"] >= age & phenoTable[,"age"] < ageMax, ]
    phenoTable <- phenoTable[phenoTable[,"race"] == race, ]

    phenoVector <- phenoTable[,dataId]
    names(phenoVector) <- phenoTable[,gtexId]
    
    return (phenoVector)

}

readGTEXmetadata <- function(x, inds = NULL) {
    # Reads the gtex covariate file given its location and selects for the individuals of interest if specified
    
    x <- read.delim(x, sep = "\t", stringsAsFactors = F, header =T)
    
    if(!is.null(inds))
        x <- x[ x[,1] %in% inds, ]
    
    return(x)
    
}
    
readPhenotype <- function(vcfFile, phenotype){

    # Helper function that given a phenoVector(whose names are individual) returns
    # a list whose first element is the same phenovector but only with enties cotained in the VCF file,
    # the second element of the returned list is a vector of the individuals
    
    if(!is.numeric(phenotype))
        stop("phenotype has to be numeric")
    # Checking sample names and selecting only the ones in the vcf file
    individuals <- names(phenotype)
    individuals <- findNames(vcfFile, individuals)
    inVcf <- checkIndividuals(vcfFile, individuals)
    if(sum(inVcf) < length(individuals)){
        warning("Not all individiuals are contained in the vcfFile")
        if(sum(inVcf) < 2)
            stop("At least two individuals are needed to be contained in the vcfFile")
    }
    names(phenotype) <- individuals
    phenotype <- phenotype[inVcf]
    individuals <- individuals[inVcf]
    
    return(list(phenotype, individuals))        
}

getPhenoVector <- function(phenotypeTableFile, mutType, age, race, gtexId, mutId, dataId) {
    
    # Returns the phenotype vector containing the mutation load values of the given mutation type(s)
    # GENOME-WIDE
    # phenotypeTableFile - string - path to mutation table
    # mutType - vector string - mutation types of interest i.e. C>T
    # gtexId - int - column in table containing the gtex id samples
    # mutId - ind - column in table containing the mutation types
    # dataId - ind - column in table containing the mutation types

    phenoTable <- getPhenoTable(phenotypeTableFile, mutType, age, race, gtexId, mutId, dataId)
    
    phenoVector <- phenoTable[,3]
    names(phenoVector) <- phenoTable[,1]
    
    return (phenoVector)
                
}

getPhenoTable <- function(phenotypeTableFile, mutType, age, race, gtexId, mutId, dataId, ageMax = 200) {
    
    # Returns the phenotype table containing the mutation load values of the given mutation type(s)
    # Columns: gtexId, mutType, value (frequency)
    # GENOME-WIDE
    # phenotypeTableFile - string - path to mutation table
    # mutType - vector string - mutation types of interest i.e. C>T
    # gtexId - int - column in table containing the gtex id samples
    # mutId - ind - column in table containing the mutation types
    # dataId - ind - column in table containing the mutation types

    phenoTable <- read.table(phenotypeTableFile, sep = "\t", stringsAsFactors = F, header = T)

    if(sum(mutType %in% phenoTable[,mutId]) !=  length(mutType))
        stop("Mutation type not found in phenotype table")

    phenoTable <- phenoTable[ phenoTable[,mutId] %in%  mutType, ]
    phenoTable <- phenoTable[phenoTable[,"age"] >= age & phenoTable[,"age"] < ageMax, ]
    phenoTable <- phenoTable[phenoTable[,"race"] %in% race, ]
    
    phenoTable <- phenoTable[, c(gtexId, mutId, dataId)]
    
    return (phenoTable)
                
}



readPhenoMatrixGenes <- function(phenoMatFile, sharedGenes, pseudo0 = 0.0001, logTransformed = T) {
    
    # Reads phenoMatrix of gene by gene mutational load. Selects only genes for which at least "sharedGenes" percentage genes
    # have a valid mutation values. For FC a valid values something not equal to one for individual tissues this is greater than the pseudo0 value
    # Remember that if the file names contains FC (for it renames look.
    
    # GENE-BY-GENE
    
    # phenoMatFile - string - path to mutation table
    # sharedGenes - real - percentage of individual having a valid value in order for a gene to be included
    # pseudo0 - pseudo0 - pesudo zero value in individual tissues
    # logTransformed - logical - true if the data is log transformed
    
    # return - data.frame - rownames genes, colnames individuals
    
    phenoM <- read.table(phenoMatFile, sep = "\t", stringsAsFactors=F, check.names = F)
    colnames(phenoM) <- gsub("\\.", "-", colnames(phenoM))

    totalInd <- ncol(phenoM)

    if (grepl("FC", basename(phenoMatFile))) {
        
        zero <- 1
        if(logTransformed)
            zero <- 0
        
        phenoM <- phenoM[ rowSums(phenoM != zero & is.finite(as.matrix(phenoM)) ) > (totalInd*sharedGenes),  ]
        cat("Fold Change table read\n")
    } else if (grepl("[Ss]un", basename(phenoMatFile), perl = T) | grepl("[Nn]on", basename(phenoMatFile), perl = T) ) {
        
        zero <- pseudo0
        #if(logTransformed)
        #   zero <- log(pseudo0)
        
        phenoM <- phenoM[ rowSums(phenoM > zero & is.finite(as.matrix(phenoM)) ) > (totalInd*sharedGenes),  ]
    cat("Individual tissue table read\n")
    } else  {
        stop("Invalid pheno matrix name, it has to contain [FC|Sun|Non]")
    }

    return(phenoM)
}


# Find Names 
findNames <- function(vcfFile,x){
    require(VariantAnnotation)
    sampleNames <- unique(samples(scanVcfHeader(vcfFile)))
    for (i in 1:length(x)) {
        pos <- grep(x[i], sampleNames)
        if(length(pos)>0){
            if (length(pos) >1){
                stop("Conflicting individual names")
            }else{
                x[i] <- sampleNames[pos]
            }
        }
    }
    
    return(x)
}

# Checks if a vector of individuals is on a vcfFile

checkIndividuals <- function (vcfFile, x){
    require(VariantAnnotation)
    return(x %in% samples(scanVcfHeader(vcfFile)))
}


# Given the genotype matrix "GL" of posterior probabilities
# this functions identifies the genotype with a posterior cutoff

getGeno <- function(genoMat, posterior = 0.9, parallel = F, ...){
    
    stopifnot(posterior > 0.5)

    #Gets info from matrix
    sizeDim <- dim(genoMat)
    nPos <- sizeDim[1]
    nInd <- sizeDim[2]
    
    #Selects positions for which all individuals have posterior that
    # passes the cutoff
    sharePos <- vector(mode = "numeric")
    genoMatLog <- genoMat >= posterior
    sharedGeno <- rowSums(genoMatLog) == nInd
    
    #Selects the matrix from those positions
    genoMat <- genoMat[sharedGeno,,]
    sizeDim <- dim(genoMat)
    nPos <- sizeDim[1]
    nInd <- sizeDim[2]
    
    # Gets the genotype for each position in each individual
    genotype <- matrix(0, nrow = nPos, ncol = nInd, dimnames= list(rownames(genoMat[,,1]), colnames(genoMat[,,1]) ) )
    #browser()
    for (i in 1:nInd){
        currentInd <- genoMat[,i,]
        #t() is used for right order of genotypes corresponding to the adequate SNPs
        genotype[,i] <-  t(currentInd) [t(currentInd) >= posterior]
    }
    
    return(genotype)
}

getGenoFromGT <- function(genoMat) {
    
    genoMatResult <- matrix(0, nrow = nrow(genoMat), ncol(genoMat), dimnames = list(rownames(genoMat), colnames(genoMat)))
    
    genoMatResult[genoMat == "1|1" | genoMat == "1/1"] <- 2
    genoMatResult[genoMat == "0|1" | genoMat == "1|0" | genoMat == "0/1" | genoMat == "1/0"] <- 1
    genoMatResult[genoMat == "."] <- NA
    
    return(genoMatResult)
}
    

#########################
# Linear model helpers
linearRegressionHelper <- function(x,y, na.rm = F, inf.rm = F, coeff = F){
    cData <- data.frame(x=x,y=y)
    if(na.rm) cData <- cData[ !is.na(cData$x) & !is.na(cData$y), ]
    if(inf.rm) cData <- cData[ !is.infinite(cData$x) & !is.infinite(cData$y), ] 
    
    
    if(sd(x) == 0 | sd(y) == 0) return (c(pval = NA, r = NA))
       
    if(coeff){
        return(c(pval = lmp(lm(x~y, data = cData)), r = cor(x,y)))
    }else{
        return(c(pval = lmp(lm(x~y, data = cData))))
    }
}

lmp <- function (modelobject, coeff = F) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    summaryLm <- summary(modelobject)
    f <- summaryLm$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    #attributes(p) <- NULL
    #if (!coeff){
        return(p)
    #}else{
       #return(c(pval = p, summaryLm$r.squared))
    #}
}

# X is ignored and only for the pruposes of lapply
linearRegressionHelperPerm <- function(x,phenotype, genotype){
    phenotype <- phenotype[sample(1:length(phenotype))]
    return (linearRegressionHelper(phenotype, genotype))
}

# Same but using matrix eQTL package

linearRegressionMatHelper <- function(x,y, na.rm = F, inf.rm = F, rsquared = T,  useModel = modelLINEAR, coovariates = SlicedData$new(), 
                      pvalCutoff = 1, errorCovariance = numeric(), 
                      outFile = tempfile(), min.pv.by.genesnp = F){
    require(MatrixEQTL)
    
    if((!is.vector(x) & !is.matrix(x)) | !is.numeric(x)) stop("x has to be a numeric vector or matrix")
    if((!is.vector(y) & !is.matrix(y)) | !is.numeric(y)) stop("y has to be a numeric vector or matrix")
    
    if(is.vector(x)) x <- t(as.matrix(x))
    if(is.vector(y)) y <- t(as.matrix(y))   
    
    if(ncol(y) != ncol(x)) stop ("Discording number of members between x and y")
    
    if(na.rm){
        x <-as.matrix(x[,colSums(is.na(x)) == 0 & colSums(is.na(y)) == 0] ) 
        y <-as.matrix(y[,colSums(is.na(x)) == 0 & colSums(is.na(y)) == 0] ) 
    }
    if(inf.rm){
        x <-as.matrix(x[,colSums(is.infinite(x)) == 0 & colSums(is.infinite(y)) == 0] ) 
        y <-as.matrix(y[,colSums(is.infinite(x)) == 0 & colSums(is.infinite(y)) == 0] ) 
    }
    
    xMat <- SlicedData$new()
    yMat <- SlicedData$new()
    
    xMat$CreateFromMatrix(x)
    yMat$CreateFromMatrix(y)
    
    results <- Matrix_eQTL_engine(snps = xMat, gene = yMat, output_file_name = outFile, cvrt = coovariates, errorCovariance = errorCovariance, 
                      pvOutputThreshold = pvalCutoff, useModel = useModel, min.pv.by.genesnp = min.pv.by.genesnp )
    
    return(results$all$eqtls)
}


removeColinear <- function(x, corCutoff = 0.8) {
    #Remove colinear variables
    
    # Recursive function that eliminates columns that are correlated
    cat("in\n")
    corX <- cor(x)

    for(i in 1:nrow(corX)) {
        if (sum(corX[i,] > corCutoff) > 1) {
            x <- x[,-i]
            x <- removeColinear(x, corCutoff)
            break
        }
    }
    return(x)
}

selectMAF <- function(x, minorAlleleFreq = 0.1) {
    
    # Given a vector of numbers retunrs TRUE if the number of lowest frequency is 
    # greated than minorAlleleFreq.
    # However returns false if there is only one number
    
    freqs <- table(floor(x + 0.5))
    
    if(length(freqs) == 1)
        return(FALSE)
    
    return( min(freqs) >= (length(x) * minorAlleleFreq))
}   

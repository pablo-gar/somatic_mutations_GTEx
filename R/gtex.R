# Helper functions to handle gtex data and annotations


library("rjson")

if(!file.exists("../config.json"))
    stop("Config file not found, are you sure you are doing 'source(x, chdir = T)'")

CONFIG <- fromJSON(file = "../config.json")

# Reads metadata from gtex
readMetadataGTEX <- function(columns = NULL, location = CONFIG$auxiliaryFiles$metadataGTEX){
    
    x <- read.delim (location, header = T, sep = "\t", stringsAsFactors = F,comment.char="#",row.names = 2)
    
    if (is.null(columns)){
        
     return(x)
    
    }else{
        
        unknownColumns <- !columns %in% colnames(x)
        if (sum(unknownColumns) > 0){
            stop(paste0("The following columns are not contained in the metadata:\n    ",
                            paste(columns[unknownColumns], collapse = " ")
                )       
            )
        }else{
            return(x[,columns, drop = F])
        }
    }
}

# Reads sample annotation from gtex
readSampleAnnotationGTEX <- function(columns = NULL, location = CONFIG$auxiliaryFiles$sampleMetadataGTEX){
    
    x <- read.delim (location, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
    
    if (is.null(columns)){
         return(x)
    }else{
        
        unknownColumns <- !columns %in% colnames(x)
        if (sum(unknownColumns) > 0){
            stop(paste0("The following columns are not contained in the metadata:\n    ",
                            paste(columns[unknownColumns], collapse = " ")
                )       
            )
        }else{
            return(x[,columns, drop = F])
        }
    }
}

# Reads PCA of genotypes
readGenotypePCA <- function(location = CONFIG$auxiliaryFiles$genotypePCA){
    
    x <- read.table (location, header = T, sep = "\t", stringsAsFactors = F, row.names = 1, check.names = F)
    return(x)
    
}

# Read all metadata
read_all_metadata <- function(sraIds, readCountFile, transDiversityFile, metadataCols =  c("SMRIN", "SMTSISCH", "SMTSPAX"), indInfo = c("AGE", "GENDER", "BMI"), sraInfo = NULL, na.rm = T) {
    
    gtexIds <- sraToGtex(sraIds, formatOut = "short")
    gtexIdsLong <- sraToGtex(sraIds, formatOut = "long")

    # Read gene counts file
    readCountTable <- read.table(readCountFile, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    
    # Read transcriptome diversity
    transDiversity <-  read.table(transDiversityFile, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

    # Read metadata
    sampleMeta <- readSampleAnnotationGTEX(metadataCols)
    indMeta <- readMetadataGTEX(indInfo)
    genotypePCA <- t(readGenotypePCA())
    
    if(!is.null(sraInfo)) {
        sra <- list()
        for(i in seq_along(sraInfo)) {
            sra_current <- querySraRunTable(sraIds, refCol = 16, queryCol = sraInfo[i])
            sra_current <- as.data.frame(sra_current)
            sra[[i]] <- sra_current
        }
        sra <- do.call(cbind, sra)
        colnames(sra) <- names(sraInfo)
    }


    # Getting ids for which we have info for all
    if(na.rm) {
        allInfo <-
            gtexIdsLong %in% rownames(sampleMeta) & 
            gtexIds %in% rownames(indMeta) & 
            gtexIds %in% rownames(genotypePCA) &
            gtexIdsLong %in% rownames(transDiversity)

        sraIds <- sraIds[allInfo]
        gtexIds<- gtexIds[allInfo]
        gtexIdsLong <- gtexIdsLong[allInfo]
    }

    # Merging tables
    results <- readCountTable[sraIds,"n_uniqueMapped", drop = F]
    results <- cbind(results, transDiversity[gtexIdsLong, , drop = F]) 
    results <- cbind(results, sampleMeta[gtexIdsLong,])
    results <- cbind(results, indMeta[gtexIds,])
    results <- cbind(results, as.data.frame(genotypePCA)[gtexIds,1:3])
    
    if (!is.null(sraInfo)) {
        results <- cbind(results, sra[sraIds,])
    }
    #results <- results[,-1]

    # Eliminate columns with all na
    if(na.rm)
        results <- results[,colSums(is.na(results)) != nrow(results)]
    
    return(results)
}

#' Gets colnames of sraid file
#' @param sraTable - path to sra Table
getSraRunTableNames <- function(sraTable = CONFIG$auxiliaryFiles$sraTable) {
    con <- file(sraTable, 'r')
    line_1 <- readLines(con, n=1)
    line_1 <- unlist(strsplit(line_1, '\t'))
    names(line_1) <- 1:length(line_1)
    close(con)
    return(line_1)
}

#' General query tool for sra ids and metadata 
#' @param x - character, vector on which queries will be applied, e.g. sra run ids
#' @param refCol - numeric, column in sraTable whose values can match to x
#' @param queryCol - numeric, column in sraTable whose values are desired
#' @param sraTable - path to sra Table
querySraRunTable <- function(x, refCol, queryCol, sraTable = CONFIG$auxiliaryFiles$sraTable, ordered = T){
    
    x_original <- x
    x <- as.character(unique(x))
    
    
    # Getting ids
    if (length(x) < 500) {
        cutFields <- paste0("cut -f ", refCol, ",", queryCol)
        pattern <- paste(x, collapse = "|")
        pattern <- paste0("'", pattern, "'")
        ids <- system(paste0(c("grep", "-E",  pattern, sraTable, "|", cutFields), collapse = " "), intern = T)
    } else {
        ids <- queryFile(x = x, refCol = refCol, queryCol = queryCol, filepath = sraTable)
    }
    ids <- read.table(textConnection(ids), sep = "\t", header = F, stringsAsFactors = F)
    
    if(refCol > queryCol)
        ids <- ids[,2:1]
    
    if(ordered) {
        
        rownames(ids) <- ids[,1]
        
        # Returning 
        result <- rep(NA, length(x_original))
        names(result) <- x_original
        present <-  names(result)[names(result) %in% rownames(ids)]
        
        
        result[present] <- ids[present, 2]
    } else {
        result <- ids
    }
        
    
    return (result)
    
}

#' Returns a vector containing the queryCol of values from x contained in refCol from file
#' It uses R instead of bash
#' @param x - character, vector on which queries will be applied, e.g. sra run ids
#' @param refCol - numeric, column in sraTable whose values can match to x
#' @param queryCol - numeric, column in sraTable whose values are desired
#' @param filepath - path to file
#' @param sep - character to use as field separator
queryFile <- function(x, refCol, queryCol, filepath, sep = "\t") {
    
    con <- file(filepath, "r")
    result <- ""
    while(length(oneline <- readLines(con, n = 1, warn = FALSE)) > 0) {
        myvector <- unlist(strsplit(oneline, sep))
        if(myvector[refCol] %in% x) {
            current_cols <- paste(myvector[sort(c(refCol, queryCol))], collapse = sep)
            result <- paste0(c(result, current_cols), collapse = "\n")
            
            x <- x[x != myvector[refCol]]
        }
        if(length(x) == 0)
            break
    }
    close(con)
    #chomps first \n
    result <- sub("\\n","", result)
    return(result)
}


#'Sra to GTEX ids
#' @param x - character vector of sra run ids
#' @param sraTable - path to sra run table
#' @param formatOut - character - format of GTEX ouptut id, short is GTEX-XXX long is GTEX-XXX-YYY-ZZZ

sraToGtex <- function(x, sraTable = CONFIG$auxiliaryFiles$sraTable, formatOut = "short") {
    
    x <- as.character(x)
    
    if(!formatOut %in% c("short", "long"))
        stop("Format has to be 'short' or 'long'")
    
    if(formatOut == "short") {
        return(querySraRunTable(x, refCol = 16, queryCol = 32, sraTable = sraTable))
    } else if(formatOut == "long") {
        return(querySraRunTable(x, refCol = 16, queryCol = 18, sraTable = sraTable))
    }
    
}

#' Sra to tissues
#' @param x - character vector of sra run ids
#' @param sraTable - path to sra run table
#' @param formatOut - character - format of GTEX ouptut id, short is GTEX-XXX long is GTEX-XXX-YYY-ZZZ

sraToTissues <- function(x, sraTable = CONFIG$auxiliaryFiles$sraTable) {
    
    x <- as.character(x)
    
    return(code_friendly_tissues(querySraRunTable(x, refCol = 16, queryCol = 21, sraTable = sraTable)))
    
}

gtexLongToTissues <- function(x, sraTable = CONFIG$auxiliaryFiles$sraTable) {
    
    x <- as.character(x)
    
    return(code_friendly_tissues(querySraRunTable(x, refCol = 18, queryCol = 21, sraTable = sraTable)))
    
}

gtexToSra <- function(x, sraTable = CONFIG$auxiliaryFiles$sraTable) {
    
    x <- as.character(x)
    
    return(querySraRunTable(x, refCol = 32, queryCol = 16, sraTable = sraTable))
    
}

gtexLongToSra <- function(x, sraTable = CONFIG$auxiliaryFiles$sraTable) {
    
    x <- as.character(x)
    
    return(querySraRunTable(x, refCol = 18, queryCol = 16, sraTable = sraTable))
    
}

tissueToGtexLong <- function(x, sraTable = CONFIG$auxiliaryFiles$sraTable) {
    
    x <- as.character(x)
    
    return(querySraRunTable(x, refCol = 21, queryCol = 18, sraTable = sraTable, ordered = F))
    
}

#Reads gtex expression data
readAllGtexExpression <- function(samples = NULL, expressionFile = CONFIG$auxiliaryFiles$expresionAllTissues) {
    
    sep <-  "\t"
    filecon <- file(expressionFile, "r")
    commentLines <- readLines(filecon, n = 2)
    
    header <- unlist(strsplit(readLines(filecon, n = 1), sep))
    
    if(is.null(samples)) {
        sampleIndex <- rep(T, length(header))
        sampleIndex[2] <- F
    } else {
        sampleIndex <- header %in% samples
    }
    
    sampleIndex[1] <- T
    
    colClasses <- rep("NULL", length(header))
    colClasses[sampleIndex] <- "numeric"
    colClasses[1] <- "character"
    
    result <- read.table(filecon, header = F, sep = sep, colClasses = colClasses) 
    colnames(result) <- header[sampleIndex]
    
    close(filecon)
    
    return(result)
    
}

#Reads gtex gene read counts
readAllGeneReadCounts <- function(samples = NULL) {
    return(readAllGtexExpression(samples, expressionFile = CONFIG$auxiliaryFiles$readsGenesAllTissues))
}

# convert long gtex ids to short ones 
gtexLongToShort <- function(x) {
     gsub('(GTEX-\\w+?)-.+', '\\1', x)
}

# Convert tissues to better names
correct_tissue_names <- function(tissues) {
    
    convert_tis <- setNames(c("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic", "Small_Intestine_Terminal_Ileum", "Muscle_Skeletal", "Breast_Mammary_Tissue", "Whole_Blood"),
                            c("Skin sun exposed", "Skin sun protected", "Small Intestine", "Muscle", "Breast", "Blood"))
    
    # renaming
    for(i in seq_along(convert_tis)) {
        tissues[tissues %in% convert_tis[i]] <- convert_tis[tissues[tissues %in% convert_tis[i]]]
    }
    
    # Eliminating under score
    tissues <- gsub("_", " ", tissues)
    
    return(tissues)
    
}

#' Convert original tissues to code-friendly tissue
code_friendly_tissues <- function(x) {
    
    x <- gsub('-', '', x)
    x <- gsub('\\(', '', x)
    x <- gsub('\\)', '', x)
    x <- gsub('\\s+', '_', x)
    
    return(x)
}


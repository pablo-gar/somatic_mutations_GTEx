library("biomaRt")

ensemblToAlias <- function(x, organism = "hsapiens_gene_ensembl"){
	ensembl <- useMart("ensembl", dataset = organism, host = "www.ensembl.org", ensemblRedirect = FALSE)
	
	output <- getBM(filters = "ensembl_gene_id", values = unique(x), attributes = c("ensembl_gene_id","hgnc_symbol"), mart = ensembl)
    output <- output[!duplicated(output$ensembl_gene_id),]
    rownames(output) <- output$ensembl_gene_id
    
	genes <- vector(mode = "character", length = length(x))
	names(genes) <- x
	genes <- output[names(genes), "hgnc_symbol"]
	
	return(genes)
}

aliasToEnsembl <- function(x, organism = "hsapiens_gene_ensembl"){
	ensembl <- useMart("ensembl", dataset = organism, host = "www.ensembl.org", ensemblRedirect = FALSE)
	
	output <- getBM(filters = "hgnc_symbol", values = unique(x), attributes = c("ensembl_gene_id","hgnc_symbol"), mart = ensembl)
    output <- output[!duplicated(output$hgnc_symbol),]
    rownames(output) <- output$hgnc_symbol
    
	genes <- vector(mode = "character", length = length(x))
	names(genes) <- x
	genes <- output[names(genes), "ensembl_gene_id"]
	
	return(genes)
}

ensemblToDescription <- function(x, organism = "hsapiens_gene_ensembl"){
	ensembl <- useMart("ensembl", dataset = organism, host = "www.ensembl.org", ensemblRedirect = FALSE)
	
	output <- getBM(filters = "ensembl_gene_id", values = unique(x), attributes = c("ensembl_gene_id", "description"), mart = ensembl)
    output <- output[!duplicated(output$ensembl_gene_id),]
    rownames(output) <- output$ensembl_gene_id
    
	genes <- vector(mode = "character", length = length(x))
	names(genes) <- x
	genes <- output[names(genes), "description"]
	
	return(genes)
}

aliasToDescription <- function(x, organism = "hsapiens_gene_ensembl"){
	ensembl <- useMart("ensembl", dataset = organism, host = "www.ensembl.org", ensemblRedirect = FALSE)
	
	output <- getBM(filters = "hgnc_symbol", values = unique(x), attributes = c("hgnc_symbol", "description"), mart = ensembl)
    output <- output[!duplicated(output$hgnc_symbol),]
    rownames(output) <- output$hgnc_symbol
    
	genes <- vector(mode = "character", length = length(x))
	names(genes) <- x
	genes <- output[names(genes), "description"]
	
	return(genes)
}

#' Returns the coordinates of a vector of genes as a data.frame
#' coordinates are 1-based
#' @param x vector with gene id/names
#' @param organism - e.g.  "hsapiens_gene_ensembl"
#' @param geneType e.g. "hgnc_symbol" or " "ensembl_gene_id"
getGeneCordinates  <- function(x, organism = "hsapiens_gene_ensembl", geneType = "hgnc_symbol") {
    
    stopifnot(is.character(x))
    x <- unique(x)
    
	ensembl <- useMart("ensembl")
	ensembl <- useDataset(organism, mart=ensembl)
	output <- getBM(filters = geneType, values = x, attributes = c(geneType, "chromosome_name", "start_position", "end_position", "strand", "description"), mart = ensembl)
    
    # only main chromosomes
    output <- output[ output$chromosome_name %in% c(1:22, "X", "Y"), ]
    output <- output[ !duplicated(output[,geneType]), ]
    
    rownames(output) <- output[, geneType]
    output <- output[x, ]
    
    rownames(output) <- x
    
    return(output)
    
}


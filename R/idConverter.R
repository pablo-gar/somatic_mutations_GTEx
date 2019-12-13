if(!file.exists("./ensembl_to_symbol.txt"))
    stop("ensembl-to-symbol conversion table not found, are you sure you are doing 'source(x, chdir = T)'")

ENSEMBL_TO_SYMBOL <- read.table("./ensembl_to_symbol.txt", sep = "\t", stringsAsFactors = F, header = T, comment.char = "")

ensemblToAlias <- function(x, organism = "hsapiens_gene_ensembl"){
	require("biomaRt")
	ensembl <- useMart("ensembl")
	ensembl <- useDataset(organism,mart=ensembl)
	
	output <- getBM(filters = "ensembl_gene_id", values = x, attributes = c("ensembl_gene_id","hgnc_symbol"), mart = ensembl)
	genes <- vector(mode = "character", length = length(x))
	names(genes) <- x
	genes[ output$ensembl_gene_id ] <- output$hgnc_symbol
	
	return(genes)
}

ensemblToAlias_offline <- function(x){
    
	output <- ENSEMBL_TO_SYMBOL[ ENSEMBL_TO_SYMBOL[,2] %in% x,]
    
	genes <- vector(mode = "character", length = length(x))
	names(genes) <- x
	genes[ output[,2] ] <-  output[,1]
	
	return(genes)
}

aliasToEnsembl <- function(x, organism = "hsapiens_gene_ensembl"){
	require("biomaRt")
	ensembl <- useMart("ensembl")
	ensembl <- useDataset(organism,mart=ensembl)
	
	output <- getBM(filters = "hgnc_symbol", values = x, attributes = c("ensembl_gene_id","hgnc_symbol"), mart = ensembl)
	genes <- vector(mode = "character", length = length(x))
	names(genes) <- x
	genes[ output$hgnc_symbol ] <- output$ensembl_gene_id
	
	return(genes)
}


library("biomaRt")
library("org.Hs.eg.db")
library("topGO")

#' performGO
#' @param gene a named numeric vector, names are gene ids and values are used for selection criteria of target genes
#' @param statistic Do _ks_ for ranked list
#' @return a list the topGOdata and fisher resuts
performGO <-  function(genes, type = "Ensembl", ontology = "BP", gene_selection_function = function(p) p < 0.05, statistic = "fisher", algorithm = "classic") {
    
    GOdata <- new("topGOdata", ontology = "BP", allGenes = genes, geneSel = gene_selection_function, annot = annFUN.org, mapping="org.Hs.eg.db", ID=type)
    test <- runTest(GOdata, algorithm = algorithm, statistic = statistic)
    
    results <- list(GOdata = GOdata, test = test)
    names(results)[2] <- statistic
    
    return(results)
    
}


#' get_GO_results
#' takes output of perfomGO and returns a data.frame with GO terms and their pvalues
get_GO_results <- function(x, test = "fisher", fdr_signif = F, bonf_signif = T, n = NULL) {
        
    if(sum(fdr_signif, bonf_signif, !is.null(n)) > 1)
        stop("Only one of these arguments can be TRUE, fdr_signif, bonf_signif, n")
        
    pvals <- score(x[[test]])
        
    if (is.null(n)) {
        result <- GenTable(x$GOdata, x[[test]], numChar = 1000) 
    } else {
        result <- GenTable(x$GOdata, x[[test]], numChar = 1000, topNodes = n) 
    }   
        
            
    result$p_bonf <- p.adjust(pvals, method = "bonferroni")[result$GO.ID]
    result$fdr <-  p.adjust(pvals, method = "fdr")[result$GO.ID]
        
    if(fdr_signif)
        result <- result[result$fdr < 0.05,]
    if(bonf_signif)
        result <- result[result$p_bonf < 0.05,]
        
    return(result)
        
}

#' Returns genes that are in each category for the top n results in the output of get_GO_results
#' @param genes vector of gene ids
#' @param x output of get_GO_results
#' @param GO_results performGO
#' @param n number of GO for which genes are desired
get_genes_GO <- function(genes, x, GO_results, n) {
    GO_id_name <- setNames(x[,"Term"], x[,"GO.ID"])
    results <- lapply(x[1:n, "GO.ID"], function(GO_id) {
                          currentGenes <- genesInTerm(GO_results[["GOdata"]], GO_id)
                          currentGenes <- currentGenes[[1]][currentGenes[[1]] %in% genes]
                          if(length(currentGenes) > 0) {
                              return(data.frame(gene = currentGenes, GO.ID = GO_id, Term = GO_id_name[GO_id]))
                          } else {
                              return(data.frame(gene = "not_found", GO.ID = GO_id, Term = GO_id_name[GO_id]))
                          }
                        })
    results <- do.call(rbind, results)
    return(results)
}

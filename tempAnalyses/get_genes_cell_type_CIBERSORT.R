# Rscript get_genes_cell_type_CIBERSORT.R NK.cells.resting TRUE 0.1 /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/TPM_Lung_allGenes.txt /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Lung_all_genes.txt /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/genes_Lung_0.1_NK.cells.resting.txt

main <- function(cmdArgs = commandArgs(T)) {
    
    cell_type <- cmdArgs[1]  #cell type name
    top <- as.logical(cmdArgs[2]) # T/F used top correlated? Or bottom correlated?
    top_percentage <- as.numeric(cmdArgs[3]) # percentage to take from the top/bottom list
    gene_expression_file <- cmdArgs[4] # the gene expression file used in cebersort
    cibersort_out <- cmdArgs[5] # the cibersort results
    out_file <- cmdArgs[6] # where to put the list of genes
    
    #gene_expression_file <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/TPM_Lung_allGenes.txt"
    #cibersort_out <- "/scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/cibersort/CIBERSORT.Output_Abs_Lung_all_genes.txt"
    #cell_type <- "NK.cells.resting"
    #top_percentage <- 0.1
    
    
    gene_expression <- read.table(gene_expression_file, sep = "\t", stringsAsFactors = F, header = T, check.names = F) 
    genes <- gene_expression[,1]
    gene_expression <- gene_expression[,-1]
    
    cibersort <- read.table(cibersort_out, sep = "\t", stringsAsFactors = F, header = T)
    gene_expression <- gene_expression[,cibersort[,1]]
    
    # Correlate genes with cell type content
    cell_type_values <- cibersort[,cell_type]
    cor_genes_cell_type <- apply(gene_expression, 1, cor, y = cell_type_values, method = "sp")
    names(cor_genes_cell_type) <- genes
    cor_genes_cell_type <- cor_genes_cell_type[!is.na(cor_genes_cell_type)]
    
    # Get top correlated
    if(top) {
        top_genes <- rev(tail(sort(cor_genes_cell_type), length(cor_genes_cell_type) * top_percentage))
    } else {
        top_genes <- rev(tail(rev(sort(cor_genes_cell_type)), length(cor_genes_cell_type) * top_percentage))
    }
    
    write.table(data.frame(gene = names(top_genes), spearman = top_genes), out_file, sep = "\t", quote = F, row.names = F, col.names = T)
    
    
}

main()

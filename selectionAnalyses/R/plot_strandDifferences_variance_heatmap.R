library("ggplot2")
library("tidyr")
library("dplyr")
library("gplots")
library("reshape")
source("../../R/ggthemes.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    out_prefix <- cmdArgs[1]
    strand_diff_files <- cmdArgs[-1]
    
    # DEBUG
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences_compiled/n6_0.0_0.7_"
    #strand_diff_files <- list.files("/scratch/users/paedugar/somaticMutationsProject/selection/strand_differences", full = T, pattern = "*raw_strand_mut_avg.txt")
    
    strand_diff_table <- read_strand_diff(strand_diff_files)
    
    # plot average number of mutations z scores
    pdf(paste0(out_prefix, "avg_heatmap.pdf"))
    do_heatmap(strand_diff_table, column = avg, order_measure = std, order_var = "C>A_TRUE", out_prefix = paste0(out_prefix, "avg_heatmap.txt"))
    dev.off()
    
    # plot sd number of mutations z scores
    pdf(paste0(out_prefix, "std_heatmap.pdf"))
    do_heatmap(strand_diff_table, column = std, order_measure = std, order_var = "C>A_TRUE", out_prefix = paste0(out_prefix, "std_heatmap.txt"))
    dev.off()
    
    # plot fcs
    pdf(paste0(out_prefix, "avg_heatmap_fc.pdf"))
    do_heatmap(strand_diff_table, column = avg, order_measure = std, order_var = "C>A_TRUE", fc = T, out_prefix = out_prefix)
    dev.off()
    
    # plot boxplots
    p <- plot_box_plots(strand_diff_table, avg)
    ggsave(paste0(out_prefix, "avg_box.pdf"), p, height = 5)
    
    p <- plot_box_plots(strand_diff_table, std)
    ggsave(paste0(out_prefix, "std_box.pdf"), p, height = 5)
    

      
}

plot_box_plots <- function(x, column) {
    
    MUTATION_COLORS <- c(`C>A`="#3BC9F3", `C>G` = "#2B2E34", `C>T` = "#FC3218", `T>A` = "#CAD0CE", `T>C` = "#9CD169", `T>G` = "#F1CAC9")
    column <- enquo(column)
    
    mat <- get_matrix(x, !!column)
    tissues <- mat$tissue
    muts <- colnames(mat)
    
    mat <- t(apply(as.matrix(mat[,-1]), 1, scale))
    mat <- as.data.frame(mat)
    colnames(mat) <- muts[-1]
                        
    mat$tissue <- tissues
    mat <- gather(mat, mut_strand, value, -tissue)
    
    ggplot(mat, aes(x = mut_strand, y = value)) +
        geom_boxplot(colour = rep(MUTATION_COLORS, each =2)) +
        #scale_colour_manual(values= rep(MUTATION_COLORS, each =2)) +
        theme_noGrid() 
}
           
           
read_strand_diff <- function(x) {
    
    results <- list()
    for (i in x) {
        tissue <- gsub("(\\w+?)-.+", "\\1", basename(i))
        current <- read.table(i, sep = "\t", stringsAsFactors = F, header = T)
        current$tissue <- tissue
        results[[i]]  <- current
    }
    
    results <- do.call(rbind,results)
    rownames(results) <- 1:nrow(results)
    return(results)
}

do_heatmap <- function(x, column, order_measure, order_var, fc = F, out_prefix) {
    
    MUTATION_COLORS <- c(`C>A`="#3BC9F3", `C>G` = "#2B2E34", `C>T` = "#FC3218", `T>A` = "#CAD0CE", `T>C` = "#9CD169", `T>G` = "#F1CAC9")
    
    #Enquoing
    column <- enquo(column)
    order_measure <- enquo(order_measure)
    
    # Getting matrix with pvalues to plot
    if(fc == TRUE) {
        mat <- get_matrix(get_fc(x, !!column), fc)
        .scale <- "none"
        .col <- colorRampPalette(c("blue", "black", "yellow"))(9)
    } else {
        mat <- get_matrix(x, !!column)
        MUTATION_COLORS <- rep(MUTATION_COLORS, each = 2)
        .scale <- "row"
        .col <- colorRampPalette(c("#FFF5F0", "#A50F15"))(20)
    }
    
    # Arriging orders
    row_orders <- get_orders(x, !!order_measure, order_var, "tissue")
    
    mat <- as.data.frame(mat)
    rownames(mat) <- mat[,1]
    mat <- as.matrix(mat[row_orders,-1])
    
    # Plot heatmap
    #gplots::heatmap.2(mat, Rowv = F, Colv = F, trace = "none", ColSideColors = MUTATION_COLORS, col = .col, scale = .scale, density.info = "none")
    
    # Save stats
    mat <- melt(mat)
    write.table(mat, out_prefix, sep = "\t", col.names = F, row.names = F, quote = F)
    
}

get_matrix <- function(x, column) {
    
    column <- enquo(column)
    
    x %>%
        mutate(mut_strand = paste0(mutation, "_", onTranscribed)) %>%
        select(mut_strand, tissue, !!column) %>%
        spread(mut_strand, !!column) 
    
}

get_orders <- function(x, order_measure, order_var, return_column) {
    
    order_measure <- enquo(order_measure)
    mat <- get_matrix(x, !!order_measure)
    mat_numeric <- as.matrix(as.data.frame(mat[,-1]))
    
    #scale by rows to order based on z scores
    mat_numeric <- t(apply(mat_numeric, 1, scale))
    rownames(mat_numeric) <- mat[,1]
    colnames(mat_numeric) <- colnames(mat)[-1]
    
    # Getting final mat order with scale values
    mat <- mat[order(-mat_numeric[,order_var]), return_column, drop = T]
    
    return(mat)
}

#' Gets the foldchange of the desire columns between transcribed and non-transcribed
get_fc <- function(x, column){
    
    column <- enquo(column)
    fc_fun <- function(x) log2(x[2] / x[1])
    
    x %>%
        group_by(mutation, tissue) %>%
        summarise(fc = fc_fun(!!column)) %>%
        mutate(onTranscribed = ".") %>%
        ungroup()
}

main()

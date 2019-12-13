# Plots a heatmap of correlations and barplots of pvalues of chromatin vs mutations
# takes output of mutation_chromatin_correlation.R

library("tidyr")
library("dplyr")
library("gplots")
library("ggplot2")
source("../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    corelationsFolder <- cmdArgs[1]
    outHeatmap <- cmdArgs[2]
    outBar <- cmdArgs[3]
    #corelationsFolder <- "/scratch/users/paedugar/somaticMutationsProject/chromatin/chromatin_mutations/tissue_plots/n6_0.0_0.7"
    
    corelationfFiles <- file.path(corelationsFolder, list.files(corelationsFolder))
    corelationfFiles <- corelationfFiles[grep(".txt", corelationfFiles)]
    
    # Read and arrange correlation values
    corelationAll <- readCorelations(corelationfFiles)
    
    ############
    # Plotting heatmap 
    plotHeatmap(corelationAll, outHeatmap)
    
    ############
    # Plotting bar plot of pvalues
    p <- plotBars(corelationAll)
    
    ggsave(outBar, p, height = 5, width = 10)
}

plotBars <- function(corelationAll) {
    
    orderBy <- "H3K9me3"
    include <- c("H3K9me3", "H3K27me3", "H3K36me3")
    barColors <- c("#ff7f7f", "#ffd87f",  "#81acef")
    chromatin_pvals <-
        corelationAll %>%
        filter(mod %in% include) %>%
        mutate(signed_pval = -log10(pval), signed_pval = ifelse(Estimate > 0, signed_pval, -signed_pval)) 
    
    tissueOrder <- chromatin_pvals[chromatin_pvals$mod == orderBy, ]
    chromatin_pvals$tissue <- factor(chromatin_pvals$tissue, ordered = T, levels = tissueOrder[order(-tissueOrder$signed_pval), "tissue"])
    chromatin_pvals$mod <- factor(chromatin_pvals$mod, ordered = T, levels = include)
    
    # Get corrected pvalue
    corrected_p <- -log10(0.05 / nrow(chromatin_pvals))
    
    p <- ggplot(chromatin_pvals, aes(x = tissue, y = signed_pval, fill = mod)) +
                geom_bar(stat = "identity", position = "dodge", colour = "black") +
                ylab("Signed -log10(p-value)") + 
                scale_fill_manual(values = barColors) +
                geom_hline(yintercept = -corrected_p, colour = "gray30", linetype = "dashed") +
                geom_hline(yintercept = corrected_p, colour = "gray30", linetype = "dashed") +
                annotate("rect", colour = NA, ymin = -corrected_p, ymax = corrected_p, xmin  = -Inf, xmax = Inf, fill = "gray50", alpha = 0.3) +
                annotate("text", label = "n.s. Bonferroni", colour = "gray30", y = 0, x = Inf, vjust = 0, hjust = 1, fontface = "bold.italic") +
                theme_grid_y() +
                coord_cartesian(ylim = c(-max(abs(chromatin_pvals$signed_pval)), max(abs(chromatin_pvals$signed_pval)))) +
                theme(legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
            
    
}

plotHeatmap <- function(corelationAll, outHeatmap) {
    
    chromatin <-
        corelationAll %>%
        select(cor, mod, tissue) %>%
        spread(mod, cor) 
    
    chromatin <- as.data.frame(chromatin)
    rownames(chromatin) <- chromatin[,1]
    chromatin <- as.matrix(chromatin[,-1])
    
    HEATMAP_GRADIENT <- colorRampPalette(c("#0066ff", "black", "#e6e600"))(19)
    heatPars <- list(x = chromatin,
                 col = HEATMAP_GRADIENT, 
                 dendrogram = "none",
                 symbreaks = T,
                 trace = "none", density = "none"
                 )

    pdf(outHeatmap)
    do.call(heatmap.2, heatPars) 
    dev.off()
}

readCorelations <- function(x) {
    
    results <- list()
    for(i in x) {
        current <- read.table(i, sep = "\t", header = T, stringsAsFactors = F)
        colnames(current)[grep("Pr", colnames(current))] <- "pval"
        current$tissue <- basename(i)
        current$tissue <- gsub("-.+", "", current$tissue)
        current <- current[current$mod != "(Intercept)",]
        
        results[[i]] <- current
    }
    
    results <- do.call(rbind, results)
    rownames(results) <- 1:nrow(results)
    
    return(results)
}

main()

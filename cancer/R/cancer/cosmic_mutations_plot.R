# Plots percentages of mutations that are also observed in the cosmic database
library("ggplot2")
library("dplyr")
source("../../../R/ggthemes.R")

main <- function(cmdArgs = commandArgs(T)) {
    
    mutationDir <- cmdArgs[1]
    maf <- cmdArgs[2]
    outPlot_correct <- cmdArgs[3]
    outPlot_incorrect <- cmdArgs[4]
    
    HEIGHT  <- 4
    WIDTH <- 12
    COLOR_LIMIT_PVAL <- 30
    
    #mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/count"
    #mutationDir <- "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/count_randomized"
    #maf <- "n6_0.0_0.7"
    #outPlot_correct <-  "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/plots/n6_0.0_0.7/percentage_overlap_all_tissues.pdf"
    #outPlot_incorrect <-  "/scratch/users/paedugar/somaticMutationsProject/cancer/cosmic_mutations/plots/n6_0.0_0.7/percentage_overlap_all_tissues_wrongMut.pdf"
    
    counts <- read_cosmic_counts(mutationDir, maf)
    counts$percent_cosmic <- counts$overlap_correct / counts$total_mut * 100
    counts$percent_cosmic_incorrect <- counts$overlap_incorrect / counts$total_mut * 10
    counts$based_covered <- counts$based_covered * 3
    counts$prob_correct_overlap <- -phyper(q = counts$overlap_correct, m = counts$total_mut, n = counts$based_covered -  counts$total_mut, k = counts$cosmic_mut, lower.tail=F, log.p=T)
    counts$prob_incorrect_overlap <- -phyper(q = counts$overlap_incorrect, m = counts$total_mut, n = counts$based_covered -  counts$total_mut, k = counts$cosmic_mut, lower.tail=F, log.p=T)
    
    # Correcting multiple hypothesis
    #counts$prob_correct_overlap <- -log10(p.adjust(10^(-counts$prob_correct_overlap), method = "BH"))
    #counts$prob_incorrect_overlap <- -log10(p.adjust(10^(-counts$prob_incorrect_overlap), method = "BH"))
    counts$prob_correct_overlap <- -log10(p.adjust(10^(-counts$prob_correct_overlap), method = "bonf"))
    counts$prob_incorrect_overlap <- -log10(p.adjust(10^(-counts$prob_incorrect_overlap), method = "bonf"))
    
    # setting limits
    counts$prob_correct_overlap [counts$prob_correct_overlap > COLOR_LIMIT_PVAL] <- COLOR_LIMIT_PVAL
    counts$prob_incorrect_overlap [counts$prob_incorrect_overlap > COLOR_LIMIT_PVAL] <- COLOR_LIMIT_PVAL
    
    # Getting number of significant samples per tissue
    p_signif <- -log10(0.05)
    n_signif_correct <- get_n_signif_samples(counts, "prob_correct_overlap", p_signif)
    n_signif_correct$percent_cosmic <- Inf
    n_signif_incorrect <- get_n_signif_samples(counts, "prob_incorrect_overlap", p_signif)
    n_signif_incorrect$percent_cosmic_incorrect <- Inf
    
    # Ordering tissues by median
    tissueOrders <- sort(tapply(counts$percent_cosmic, counts$tissue, median))
    counts$tissue <- factor(counts$tissue, ordered = T, levels = names(tissueOrders))
    
    cols <-  colorRampPalette(c("#666666","#cec825","#ceb425","#ce9825","#ce7325","#ce2525"))(20)
    
    p <- ggplot(counts, aes(x = tissue, y = percent_cosmic)) +
            geom_jitter(aes(colour = prob_correct_overlap), width = 0.3, alpha = 0.6)+
            geom_boxplot(lwd = 0.6, width = 0.3, colour = "#00c0c6", fill = "grey80", alpha = 0) +
            geom_text(aes(label = label), data = n_signif_correct, vjust = 1, size = 1.8, angle = 45) +
            xlab("Tissue") +
            ylab("Mutations seen in cancer(%)") +
            labs(colour = "-log10(Bonferroni p-value)") +
            #coord_cartesian(ylim = c(min(counts$percent_cosmic), max(counts$percent_cosmic))) +
            coord_cartesian(ylim = c(0, 25)) +
            #scale_colour_gradientn(colours = cols, limits = c(0,max(counts$prob_correct_overlap))) +
            scale_colour_gradientn(colours = cols, limits = c(0,COLOR_LIMIT_PVAL)) +
            theme_grid_y() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1))
        #p
        
    pincorrect <- ggplot(counts, aes(x = tissue, y = percent_cosmic_incorrect)) +
            geom_jitter(aes(colour = prob_incorrect_overlap), width = 0.3, alpha = 0.6)+
            geom_boxplot(lwd = 0.6, width = 0.3, colour = "#00c0c6", fill = "grey80", alpha = 0) +
            geom_text(aes(label = label), data = n_signif_incorrect, vjust = 1, size = 1.8, angle = 45) +
            xlab("Tissue") +
            ylab("Position seen in cancer, but different mutation(%)") +
            labs(colour = "-log10(Bonferroni p-value)") +
            #coord_cartesian(ylim = c(min(counts$percent_cosmic), max(counts$percent_cosmic))) +
            coord_cartesian(ylim = c(0, 25)) +
            #scale_colour_gradientn(colours = cols, limits = c(0,max(counts$prob_correct_overlap))) +
            scale_colour_gradientn(colours = cols, limits = c(0,COLOR_LIMIT_PVAL)) +
            theme_grid_y() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1))
        #pincorrect
    
    ggsave(outPlot_correct, p, height = HEIGHT, width = WIDTH)
    ggsave(outPlot_incorrect, pincorrect, height = HEIGHT, width = WIDTH)
    
    
    
}

read_cosmic_counts <- function (mutationDir, maf) {
    
    results <- list()
    i <- 1 
    for(tissue in list.dirs(mutationDir, recursive = F)) {
        if (grepl("EXO", tissue))
            next
        for(current_sample in list.files(file.path(tissue, maf))) {
            current <- read.table(file.path(tissue, maf, current_sample), sep = "\t", stringsAsFactors = F, header = F)
            colnames(current) <- c("ind", "overlap_correct", "overlap_incorrect", "total_mut", "cosmic_mut", "based_covered")
            current$tissue <- basename(tissue)
            results[[i]] <- current
            i <- i + 1
        }
    }
    
    results <- do.call(rbind, results)
    
    return(results)
    
}

get_n_signif_samples <- function(counts, column, p_signif) {
    
    n_signif <- data.frame(total = tapply(counts[,column], counts$tissue, length),
                                   significant = tapply(counts[,column], counts$tissue, function(x) sum(x > p_signif))
                                   )
    n_signif$tissue <- rownames(n_signif)
    n_signif$percentage <- round(100 * n_signif$significant / n_signif$total)
    n_signif$label <- paste0(n_signif$percentage, "%\n(", n_signif$significant, "/", n_signif$total, ")")
    
    return(n_signif)
    
}
    

main()

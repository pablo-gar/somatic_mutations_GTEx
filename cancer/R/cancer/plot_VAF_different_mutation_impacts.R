library("ggplot2")
library("dplyr")
source("../../../R/ggthemes.R", chdir = T)
source("../../../R/gtex.R", chdir = T)
source("../../../R/misc.R", chdir = T)
source("../../../R/mutationGeneAnnotation.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    unique_mutations <- as.logical(cmdArgs[1])
    plot_last_exon <- as.logical(cmdArgs[2])
    out_prefix <- cmdArgs[3]
    mutations_file <- cmdArgs[4]
    
    # INCLUDING LAST EXON
    #mutations_file <- "/scratch/users/paedugar/somaticMutationsProject/cancer/dndsout_tissuesMerged_all/1_all_mutations_annotated_dnds.txt.gz"
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/cancer/VAF_missense_nonsense/median_vaf_per_impact_last_exon"
    #plot_last_exon <- T
    # NOT LAST EXON
    #out_prefix <- "/scratch/users/paedugar/somaticMutationsProject/cancer/VAF_missense_nonsense/median_vaf_per_impact"
    #plot_last_exon <- F
    
    
    dir.create(dirname(out_prefix), recursive = T)
    
    mutations <- read.table(mutations_file, sep = "\t", stringsAsFactors = F, header = T)
    mutations$VAF <- mutations$n_support / mutations$coverage
    
    ## Get only singletons
    #mutations <- get_mutations_n(mutations, n_muts = 1)
    
    # Annotate if a mutation is the last exon
    if(plot_last_exon) {
        mutations <- annotate_last_exon(mutations)
        impact_types <- c("Synonymous", "Missense", "NonsenseTRUE", "NonsenseFALSE")
        mutations[mutations$impact == "Nonsense", "impact"] <- paste0(mutations[mutations$impact == "Nonsense", "impact"], mutations[mutations$impact == "Nonsense", "last_exon"])
    } else {
        impact_types <-  c("Synonymous", "Missense", "Nonsense")
    }
    
    
    mutations <- mutations[ mutations$impact %in% c("Synonymous", "Missense", "Nonsense"), ]
    mutations$impact <- factor(mutations$impact, levels = impact_types, ordered = T)
    
    
    # Get unique mutations, for those repeated get their median VAFs
    if(unique_mutations) {
        mutations <- make_unique(mutations)
    }
    
    
    # ANALYSIS OF VAF PER POSITION IN EXON
    # Append where the mutation occur relative to the total size of the gene
    #mutations <- annotate_position(mutations)
    #
    #
    #bbb <- mutations[grep("Mis", mutations$impact),]
    ##bbb <- bbb[!bbb$last_exon,]
    #bbb <- tapply(bbb$VAF, bbb$gene_pos_bracket, median)
    #toPlot <- data.frame(gene_pos = names(bbb), median_VAF = bbb)
    #p <- ggplot(toPlot, aes(x = gene_pos, y = median_VAF)) +
    #    geom_bar(stat = "identity") +
    #    theme_grid_y() 
    #
    #ggsave(paste0(out_prefix, "_nonSense_vaf_across_gene.pdf"), p)
    
    
    
    # Append tissue
    samples <- unique(mutations$sample)
    samples <- setNames(sraToTissues(samples), samples)
    mutations$tissue <- samples[mutations$sample]
    brain_tisssues <- grepl("Brain",  mutations$tissue)
    
    
    
    mutations_all <- mutations
    # Getting median values and confidence intervals
    for(tis in c("All", "Brain_sections", "Non-brain")) {
        
        if(tis == "Brain_sections") {
            mutations <- mutations_all[brain_tisssues,]
        } else if (tis == "All") {
            mutations <- mutations_all
        } else {
            mutations <- mutations_all[!brain_tisssues,]
            
        }
        
        VAF <- list()
        impacts <- unique(mutations$impact)
            
        for(i in seq_along(impacts)) {
            current <- mutations[mutations$impact == impacts[i],]
            conf_int <- bootstrap_confidence_interval(current$VAF, FUN = median, bootstrap_counts = 1000)
            result <- data.frame(VAF_med = median(current$VAF), low_conf = conf_int[1], high_conf = conf_int[2], impact = impacts[i], stringsAsFactors = F)
            VAF[[i]] <- result
        }
        
        VAF <- do.call(rbind, VAF)
        VAF$impact <- factor(VAF$impact, levels = impact_types, ordered = T)
        
        # Getting wilcox pvalues
        wilcox <- list()
        impacts <- unique(mutations$impact)
        counter <- 1
        for(i in seq_along(impacts)[-length(impacts)]) {
            for(j in (i+1):length(impacts)) {
                impact_a <- impacts[i]
                impact_b <- impacts[j]
                wilcox_results <- wilcox.test(mutations[mutations$impact == impact_a, "VAF"], mutations[mutations$impact == impact_b, "VAF"])
                wilcox[[counter]] = data.frame(impact_a = impact_a, impact_b = impact_b, pvalue = wilcox_results[["p.value"]], stringsAsFactors = F)
                counter <- counter + 1
            }
        }
        wilcox <- do.call(rbind, wilcox)
        wilcox$pvalue_str <- ifelse(wilcox$pvalue == 0, "p < 1e-300", paste0("p = ", signif(wilcox$pvalue,2)))
        
        if (plot_last_exon) {
            wilcox$x <- c(1.5, 3, 2.5, 2.5, 2, 3.5)
            wilcox$y <- c(VAF[as.character(VAF$impact) == "Synonymous", "VAF_med"], 
                          VAF[as.character(VAF$impact) == "Missense", "VAF_med"], 
                          VAF[as.character(VAF$impact) == "Missense", "VAF_med"], 
                          VAF[as.character(VAF$impact) == "Synonymous", "VAF_med"],
                          VAF[as.character(VAF$impact) == "Synonymous", "VAF_med"],
                          VAF[as.character(VAF$impact) == "NonsenseTRUE", "VAF_med"])
            
            barCols <- c("#97BC7F", "#DB7171", "#8459C5", "#a992cc")
            
            ylim <-  c(0, 0.12)
        } else {
            wilcox$x <- c(1.5, 2, 2.5)
            wilcox$y <- c(VAF[as.character(VAF$impact) == "Synonymous", "VAF_med"], 
                          VAF[as.character(VAF$impact) == "Synonymous", "VAF_med"],
                          VAF[as.character(VAF$impact) == "Missense", "VAF_med"]
                         ) 
            barCols <- c("#97BC7F", "#DB7171", "#8459C5")
            ylim <-  c(0, 0.08)
        }
        
        p_median <- ggplot(VAF, aes(x = impact, y = VAF_med)) +
            geom_bar(stat = "identity", fill = barCols, colour = "black") +
            geom_errorbar(aes(ymax = high_conf, ymin = low_conf), width = 0.3) +
            geom_text(aes(x = x, y = y, label = pvalue_str), data = wilcox, vjust = 0, fontface = "italic") + 
            ylab("Median VAF") + 
            coord_cartesian (ylim = ylim) +
            #geom_text(aes(label = pvalue), fontface = "italic", vjust = 0) + 
            theme_grid_y() + 
            theme(axis.text.x = element_text(angle = 20, hjust = 1))
        
        ggsave(paste0(out_prefix, tis, ".pdf"), p_median, width = 2.5, height = 5)
        
        #p_box <- ggplot(mutations, aes(x = impact, y = VAF)) +
        #    geom_boxplot(fill = c("#97BC7F", "#DB7171", "#8459C5"), colour = "black") +
        #    ylab("Median VAF") + 
        #    theme_grid_y()
    }
        
        
}


annotate_last_exon <- function(mutations) {
    
    mutations$last_exon <- F
    last_exon <- get_last_exon_per_genes()
    rownames(last_exon) <- last_exon$gene_name
    
    last_exon <- last_exon[ last_exon$gene_name %in% mutations$gene_name, ]
    
    for(gene_name in last_exon$gene_name) {
        
        start_exon <- last_exon[gene_name, "last_exon_start"]
        end_exon <- last_exon[gene_name, "last_exon_end"]
        
        pos_to_modify <- mutations$gene_name == gene_name & mutations$pos >= start_exon & mutations$pos <= end_exon
        
        if(any(pos_to_modify))
            mutations[ pos_to_modify, "last_exon"] <- T
        
    }
    
    return(mutations)
    
}

annotate_position <- function(mutations, bins = 7) {
    
    mut_pos <- as.numeric(gsub("[A-Z](\\d+)[A-Z]", "\\1", toupper(mutations$ntchange)))
    mutations$gene_pos <- mut_pos/mutations$cds_length
    
    bin_limits <- seq(0,1, length.out = bins)
    label <- paste0(signif(bin_limits[-bins], 2), "_", signif(bin_limits[-1],2))
    
    mutations$gene_pos_bracket <- cut(mutations$gene_pos, breaks = bin_limits, labels = label, include.lowest = T)
    
    return(mutations)
}

get_mutations_n <- function(mutations, n_muts) {
    
    ids <- with(mutations, paste0(chr, ".", pos, ".", ref, mut))
    
    id_counts <- table(ids)
    
    mutations <- mutations[ ids %in% names(id_counts)[id_counts == n_muts], ]
    
    return(mutations)
    
}
    
make_unique <- function(mutations) {
    
    id <- with(mutations, paste0(chr, ".", pos, ".", ref, mut, ".", gene_name, ".", impact))
    
    medians <- tapply(mutations$VAF, id, median)
    mutations <- mutations[!duplicated(id),]
    
    id <- with(mutations, paste0(chr, ".", pos, ".", ref, mut, ".", gene_name, ".", impact))
    mutations$VAF <- as.vector(medians[id])
    
    return(mutations)
    
}

main()

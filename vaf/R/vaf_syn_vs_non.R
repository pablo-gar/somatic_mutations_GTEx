# Performs basic analysis of VAF and metadata
source("../../R/plots.R", chdir=T)
source("../../R/ggthemes.R", chdir=T)
library("dplyr")
library("tidyr")
library("ggplot2")
library('ggrepel')

main <- function(cmdArgs = commandArgs(T)) {
    
    mutation_map_file <- cmdArgs[1]
    out_folder <- cmdArgs[2]
    mutation_map_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/all_mutations_annotated.txt"
    metadata_file <- "/scratch/users/paedugar/somaticMutationsProject/mutationCount/1_final_tables/metadata.txt"
    out_folder <- "/scratch/users/paedugar/somaticMutationsProject/vaf/syn_vs_non_vaf_analysis/"
    
    ## READ mutations
    mutation_map <- read.table(mutation_map_file, sep="\t", stringsAsFactors=F, header=T)
    
    ## READ metadata
    metadata <- read.table(file=metadata_file, header=T, stringsAsFactors=F)
    
    vaf_per_sample <- get_vaf_per_sample_impact(mutation_map, vaf, metadata)
    #vaf_per_sample_corrected <- get_vaf_per_sample_impact(mutation_map, vaf, metadata)
    
    #---------------------------------------# 
    # Basic correaltions
    #---------------------------------------# 
    
    #syn/mis vs syn/non
    p <- scatter(as.data.frame(vaf_per_sample), x='syn_mis', y='syn_non', pSize=1.5, alpha=0.3) +
        theme_sleek() +
        xlab("log2(Synonymous VAF / Missense VAF)") + ylab("log2(Synonymous VAF / Nonsense VAF)") + 
        theme_grid()
    
    ggsave(paste0(out_folder, "syn_mis_vs_syn_non.pdf"), p)
    
    #syn/non vs seq depth
    p <- scatter(as.data.frame(vaf_per_sample), y='syn_non', x='n_uniqueMapped', pSize=1.5, alpha=0.3) +
        ylab("log2(Synonymous VAF / Nonsense VAF)") + xlab("Sequencing depth") + 
        theme_sleek() +
        theme_grid()
    
    ggsave(paste0(out_folder, "syn_non_vs_seq_depth.pdf"), p)
    
    #syn/mis vs seq depth
    p <- scatter(as.data.frame(vaf_per_sample), y='syn_mis', x='n_uniqueMapped', pSize=1.5, alpha=0.3) +
        ylab("log2(Synonymous VAF / Missense VAF)") + xlab("Sequencing depth") + 
        theme_sleek() +
        theme_grid()
    
    ggsave(paste0(out_folder, "syn_mis_vs_seq_depth.pdf"), p)
    
    # sys/non vs transcriptome diversity
    p <- scatter(as.data.frame(vaf_per_sample), y='syn_non', x='transcriptome_diversity', pSize=1.5, alpha=0.3) +
        ylab("log2(Synonymous VAF / Nonsense VAF)") + xlab("Transcriptome diversity") +
        theme_sleek()
    
    ggsave(paste0(out_folder, "syn_non_vs_transcriptome_diversity.pdf"), p)
    
    # all metadata cors
    p_syn_mis_metadata_cors <- plot_cors_metadata(vaf_per_sample, syn_mis) 
    p_syn_non_metadata_cors <- plot_cors_metadata(vaf_per_sample, syn_non)
    
    ggsave(paste0(out_folder, "metadata_cors_syn_mis.pdf"), p_syn_mis_metadata_cors, width=14)
    ggsave(paste0(out_folder, "metadata_cors_syn_non.pdf"), p_syn_non_metadata_cors, width=14)
    
    #---------------------------------------# 
    # VAF analysis for synonymous vs non synonymous
    
    # Get median VAF per sample along tissues
    
    p <- plot_metric_tissues(vaf_per_sample, syn_mis)
    ggsave(paste0(out_folder, "syn_mis_tissues.pdf"), p, width=14)
    
    p <- plot_metric_tissues(vaf_per_sample, syn_non)
    ggsave(paste0(out_folder, "syn_non_tissues.pdf"), p, width=14)
    
    # Plot those medians vs other metadata
    p <- plot_metric_tissues_seqDepth(vaf_per_sample, syn_mis)
    ggsave(paste0(out_folder, "syn_mis_tissues_seq_depth.pdf"), p)
    
    p <- plot_metric_tissues_seqDepth(vaf_per_sample, syn_non)
    ggsave(paste0(out_folder, "syn_non_tissues_seq_depth.pdf"), p)
    #---------------------------------------# 
    
}

get_vaf_per_sample_impact <- function(x, metric, metadata) {
    
    metric <- enquo(metric)
    
    vaf_per_sample <- x %>%
        filter(!is.na(impact)) %>%
        group_by(sraIds, impact) %>%
        summarise(gtexIds = gtexIds[1], tissue=tissue[1], vaf_median = median(!!metric)) %>%
        ungroup() %>%
        spread(impact, vaf_median) %>%
        mutate(syn_mis = log2(Synonymous/Missense), syn_non = log2(Synonymous/Nonsense)) %>%
        filter(!is.na(syn_mis)) %>%
        mutate(tissue = factor(tissue, ordered=T, levels=names(sort(tapply(syn_mis, tissue, median)))))
    
    vaf_per_sample <- left_join(vaf_per_sample, metadata, by="sraIds")
    
    return(vaf_per_sample)
}

plot_cors_metadata <- function(x, metric) {
    
    metric <- enquo(metric)
        
    toPlot <- gather(x, "feature", "value", n_uniqueMapped, transcriptome_diversity, AGE, GENDER, BMI, C1, C2)
    spearman_cors <- toPlot %>%
    group_by(tissue, feature) %>%
    summarise(spearman=cor.test(value, !!metric, method="spearman")[["estimate"]], 
              pvalue=cor.test(value, !!metric, method="spearman")[["p.value"]],
              log10_signed_pvalue = -log10(pvalue) *(abs(spearman)/ spearman)) %>%
    ungroup() %>%
    mutate(tissue=factor(tissue, ordered=T, levels=names(sort(tapply(log10_signed_pvalue[feature=="AGE"], tissue[feature=="AGE"], median)))))
    
    ggplot(spearman_cors, aes(x=log10_signed_pvalue, y=tissue)) +
        geom_point() +
        geom_vline(xintercept=0, colour="grey50", linetype='dashed') + 
        facet_grid(.~feature) +
        theme_bw() + 
        theme_sleek() + 
        theme_grid()
        
}

plot_metric_tissues <- function(x, metric) {
    
    metric <- enquo(metric)
    toPlot <- x %>% 
        dplyr::filter(!is.na(!!metric)) %>%
        mutate(tissue = factor(tissue, ordered=T, levels=names(sort(tapply(!!metric, tissue, median)))))
    
    
    ggplot(toPlot, aes(x = tissue, y = !!metric)) +
        geom_boxplot(width=0.3) +
        geom_violin(fill=NA) +
        theme_sleek() +
        theme_grid_y() +
        theme(axis.text.x = element_text(angle=30, hjust=1))
}


plot_metric_tissues_seqDepth <- function(x, metric) {
    
    metric_or <- deparse(substitute(metric))
    metric <- enquo(metric)
    
    toPlot <- x %>% 
        dplyr::filter(!is.na(!!metric) & !is.na(n_uniqueMapped)) %>%
        group_by(tissue) %>%
        summarise(value=median(!!metric), seq_depth=sum(as.numeric(n_uniqueMapped))) %>%
        ungroup()
    
    scatter(as.data.frame(toPlot), y="value", x="seq_depth", alpha=1) + 
        geom_text_repel(aes(label=tissue)) +
        xlab('Sequencing depth') + ylab(metric_or) +
        theme_sleek()
}

main()
